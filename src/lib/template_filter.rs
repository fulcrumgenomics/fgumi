//! The pre-grouping template filter shared by the `group` and `dedup` commands.
//!
//! Both commands drop templates that fail the same criteria before grouping them
//! by position and UMI. This module owns the single implementation; per-reason
//! accounting goes into a [`TemplateFilterCounts`], which each command renders
//! into its own metrics schema (see `fgumi_metrics::template_filter`).
//!
//! The two commands were previously separate, near-identical copies that differed
//! only in how they handled a fully unmapped template — the behaviour now
//! expressed by [`TemplateFilterConfig::allow_unmapped`].

use crate::metrics::{TemplateFilterCounts, TemplateFilterReason};
use crate::sam::SamTag;
use crate::template::Template;
use crate::umi::{UmiValidation, validate_umi};
use fgumi_raw_bam;
use fgumi_raw_bam::RawRecord;

/// Criteria for [`filter_template`].
#[derive(Debug, Clone)]
pub struct TemplateFilterConfig {
    /// UMI tag bytes (e.g. `[b'R', b'X']`).
    pub umi_tag: [u8; 2],
    /// Minimum mapping quality, applied to a read's own `MAPQ` and to its mate's
    /// `MQ` tag.
    pub min_mapq: u8,
    /// Keep reads flagged QC-fail rather than discarding them.
    pub include_non_pf: bool,
    /// Minimum UMI length; `None` disables the check.
    pub min_umi_length: Option<usize>,
    /// Skip UMI validation entirely (missing UMI, Ns, and length are not checked).
    pub no_umi: bool,
    /// Accept templates whose every primary read is unmapped.
    ///
    /// `group` sets this from `--allow-unmapped`. `dedup` always passes `false`:
    /// its `--include-unmapped` pass-throughs are split off *before* the filter
    /// runs, so an unmapped template reaching this function is always a rejection.
    pub allow_unmapped: bool,
}

/// Returns `true` when any record in the template is shorter than the minimum BAM
/// record length, i.e. the input is truncated or corrupt.
///
/// Shared with `dedup`'s `--include-unmapped` pass-through check so the two agree
/// on what "malformed" means: a corrupt record must never be passed through.
#[must_use]
pub fn template_has_malformed_record(template: &Template) -> bool {
    template.records().iter().any(|r| r.len() < fgumi_raw_bam::MIN_BAM_RECORD_LEN)
}

/// Returns `true` when the template has no mapped primary read.
///
/// Shared with `dedup`'s `--include-unmapped` pass-through check so the two agree
/// on what "unmapped" means. A template with neither primary read present is not
/// "unmapped" — see [`TemplateFilterReason::NoPrimaryReads`] — so this returns
/// `false` for it.
#[must_use]
pub fn template_is_fully_unmapped(template: &Template) -> bool {
    let (raw_r1, raw_r2) = (template.r1(), template.r2());
    if raw_r1.is_none() && raw_r2.is_none() {
        return false;
    }
    raw_r1.is_none_or(RawRecord::is_unmapped) && raw_r2.is_none_or(RawRecord::is_unmapped)
}

impl Default for TemplateFilterConfig {
    /// The CLI defaults: standard `RX` UMI tag, no mapping-quality or UMI-length
    /// minimum, and no exceptions. Deliberately hand-written rather than derived —
    /// a derived `umi_tag` would be `[0, 0]`, which is not a valid tag.
    fn default() -> Self {
        Self {
            umi_tag: *SamTag::RX,
            min_mapq: 0,
            include_non_pf: false,
            min_umi_length: None,
            no_umi: false,
            allow_unmapped: false,
        }
    }
}

/// Filters a template against `config`, recording the outcome in `counts`.
///
/// Returns `true` when the template should be kept. Every call records exactly
/// one outcome — one accept or one rejection with a specific
/// [`TemplateFilterReason`] — so the counts always reconcile.
pub fn filter_template(
    template: &Template,
    config: &TemplateFilterConfig,
    counts: &mut TemplateFilterCounts,
) -> bool {
    // Primary-record count, recorded alongside the template count so `group` and
    // `dedup` can render the same measurement in their respective units.
    let primary_reads = u64::from(template.r1().is_some()) + u64::from(template.r2().is_some());

    // Fail closed when *any* record (primary, secondary, or supplementary) is
    // shorter than the minimum BAM record length: a truncated record indicates
    // corrupt input, so the template must be rejected rather than have the
    // malformed record silently dropped (and later panic in RawRecordView::new
    // on the dedup/serialize path).
    if template_has_malformed_record(template) {
        counts.record_rejected(TemplateFilterReason::MalformedRecord, primary_reads);
        return false;
    }

    let raw_r1 = template.r1();
    let raw_r2 = template.r2();

    if raw_r1.is_none() && raw_r2.is_none() {
        counts.record_rejected(TemplateFilterReason::NoPrimaryReads, primary_reads);
        return false;
    }

    if template_is_fully_unmapped(template) && !config.allow_unmapped {
        counts.record_rejected(TemplateFilterReason::Unmapped, primary_reads);
        return false;
    }

    // Phase 1: Flag-based checks
    for raw in [raw_r1, raw_r2].into_iter().flatten() {
        if !config.include_non_pf && raw.is_qc_fail() {
            counts.record_rejected(TemplateFilterReason::NotPassingFilter, primary_reads);
            return false;
        }

        if !raw.is_unmapped() {
            let mapq = fgumi_raw_bam::mapq(raw);
            if mapq < config.min_mapq {
                counts.record_rejected(TemplateFilterReason::LowMappingQuality, primary_reads);
                return false;
            }
        }
    }

    // Phase 2: Single-pass tag lookups (MQ + UMI in one aux scan)
    for raw in [raw_r1, raw_r2].into_iter().flatten() {
        let aux = fgumi_raw_bam::aux_data_slice(raw);
        let check_mq = !raw.is_mate_unmapped();
        let check_umi = !config.no_umi;

        let (found_mq, found_umi) =
            scan_aux_for_mq_and_umi(aux, config.umi_tag, check_mq, check_umi);

        // Compare as signed so a negative MQ (e.g. `MQ:c:-1`) fails the filter rather
        // than wrapping to 255 via `as u8`. `found_mq` is `Some` only when `check_mq`
        // was set, so no further guard is needed here.
        if let Some(mq) = found_mq
            && mq < i64::from(config.min_mapq)
        {
            counts.record_rejected(TemplateFilterReason::LowMateMappingQuality, primary_reads);
            return false;
        }

        // Skip UMI validation in no-umi mode
        if config.no_umi {
            continue;
        }

        if let Some(umi_bytes) = found_umi {
            match validate_umi(umi_bytes) {
                UmiValidation::ContainsN => {
                    counts.record_rejected(TemplateFilterReason::NsInUmi, primary_reads);
                    return false;
                }
                UmiValidation::Valid(base_count) => {
                    if let Some(min_len) = config.min_umi_length
                        && base_count < min_len
                    {
                        counts.record_rejected(TemplateFilterReason::UmiTooShort, primary_reads);
                        return false;
                    }
                }
            }
        } else {
            counts.record_rejected(TemplateFilterReason::MissingUmi, primary_reads);
            return false;
        }
    }

    counts.record_accepted(primary_reads);
    true
}

/// Scans a record's aux block once for both the `MQ` tag and the UMI tag.
///
/// Returns `(mate_mapping_quality, umi_bytes)`, either of which is `None` when the
/// tag is absent, unparseable, or was not requested. A single pass is used because
/// this runs once per primary read of every template in the input.
///
/// `check_mq` and `check_umi` let the caller skip the half it does not need — the
/// scan stops as soon as everything requested has been found.
fn scan_aux_for_mq_and_umi(
    aux: &[u8],
    umi_tag: [u8; 2],
    check_mq: bool,
    check_umi: bool,
) -> (Option<i64>, Option<&[u8]>) {
    if !check_mq && !check_umi {
        return (None, None);
    }
    let mut found_mq: Option<i64> = None;
    let mut found_umi: Option<&[u8]> = None;
    let mut p = 0;
    while p + 3 <= aux.len() {
        let t = [aux[p], aux[p + 1]];
        let val_type = aux[p + 2];

        if check_umi && t == umi_tag && val_type == b'Z' {
            let start = p + 3;
            if let Some(end) = aux[start..].iter().position(|&b| b == 0) {
                found_umi = Some(&aux[start..start + end]);
                p = start + end + 1;
            } else {
                break;
            }
            if !check_mq || found_mq.is_some() {
                break;
            }
            continue;
        }

        if check_mq && t == *SamTag::MQ {
            found_mq = fgumi_raw_bam::extract_int_value(aux, p, val_type);
        }

        if let Some(size) = fgumi_raw_bam::tag_value_size(val_type, &aux[p + 3..]) {
            p += 3 + size;
        } else {
            break;
        }
        if (!check_umi || found_umi.is_some()) && (!check_mq || found_mq.is_some()) {
            break;
        }
    }
    (found_mq, found_umi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::SamBuilder;

    /// Builds a template whose every primary read is unmapped.
    fn unmapped_template() -> Template {
        let mut builder = SamBuilder::new();
        builder
            .read_name(b"unmapped")
            .sequence(b"ACGT")
            .qualities(&[30; 4])
            .flags(fgumi_raw_bam::flags::UNMAPPED | fgumi_raw_bam::flags::MATE_UNMAPPED)
            .add_string_tag(SamTag::RX, b"ACGT");
        Template::from_records(vec![builder.build()]).expect("template construction")
    }

    fn config(allow_unmapped: bool) -> TemplateFilterConfig {
        TemplateFilterConfig {
            umi_tag: *SamTag::RX,
            min_mapq: 0,
            include_non_pf: false,
            min_umi_length: None,
            no_umi: false,
            allow_unmapped,
        }
    }

    /// `allow_unmapped` is the only behavioural difference between the two callers
    /// this function was unified from: `group` sets it from `--allow-unmapped`,
    /// while `dedup` always passes `false` because its `--include-unmapped`
    /// pass-throughs are split off before the filter runs. Both settings are
    /// exercised here so the merge cannot silently lose one command's behaviour.
    #[rstest::rstest]
    #[case::rejected_when_disallowed(false, false)]
    #[case::accepted_when_allowed(true, true)]
    fn unmapped_templates_respect_allow_unmapped(
        #[case] allow_unmapped: bool,
        #[case] expect_accept: bool,
    ) {
        let mut counts = TemplateFilterCounts::new();
        let accepted = filter_template(&unmapped_template(), &config(allow_unmapped), &mut counts);

        assert_eq!(accepted, expect_accept);
        assert_eq!(
            counts.rejected_templates(TemplateFilterReason::Unmapped),
            u64::from(!expect_accept)
        );
        assert_eq!(counts.accepted_templates(), u64::from(expect_accept));
    }
}
