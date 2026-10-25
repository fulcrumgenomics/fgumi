//! The alignment scoring an aligner command line aligns with.
//!
//! The realign stage recomputes `AS` for an alignment it promotes from `XA`, so it
//! needs the scoring the aligner used: [`BwaScoring::from_command`] reads it from a
//! subprocess command line (a preset's or `--aligner::command`'s), and the
//! in-process backend reads it from its own options.

/// bwa-mem3 alignment scoring, used to recompute `AS` for a promoted alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BwaScoring {
    /// Score for a matching base (`-A`).
    pub(crate) a: i32,
    /// Penalty for a mismatching base (`-B`).
    pub(crate) b: i32,
    /// Deletion open penalty (`-O`, first value).
    pub(crate) o_del: i32,
    /// Deletion extension penalty (`-E`, first value).
    pub(crate) e_del: i32,
    /// Insertion open penalty (`-O`, second value).
    pub(crate) o_ins: i32,
    /// Insertion extension penalty (`-E`, second value).
    pub(crate) e_ins: i32,
}

impl Default for BwaScoring {
    /// bwa-mem3 defaults: `-A 1 -B 4 -O 6,6 -E 1,1`.
    fn default() -> Self {
        Self { a: 1, b: 4, o_del: 6, e_del: 1, o_ins: 6, e_ins: 1 }
    }
}

/// bwa-mem3's `--meth-scoring` modes. Only `Collapsed` changes the scoring this
/// module models (its `-B 2`); `Genomic` and `Neutral` keep bwa's `-B`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MethScoring {
    Collapsed,
    Genomic,
    Neutral,
}

/// What a command's words set, before defaults are applied.
#[derive(Debug, Default)]
struct ScoringWords<'a> {
    a: Option<i32>,
    b: Option<i32>,
    open: Option<(i32, i32)>,
    extend: Option<(i32, i32)>,
    /// `-x <mode>`.
    mode: Option<&'a str>,
    /// `--meth` (`Some(false)`) or `--meth=taps` (`Some(true)`).
    meth_taps: Option<bool>,
    meth_scoring: Option<MethScoring>,
}

impl BwaScoring {
    /// The mismatch penalty `bwa-mem3 mem --meth` sets under its default
    /// collapsed bisulfite scoring (bwameth's `-B 2`), scaled by `-A`.
    const METH_MISMATCH: i32 = 2;

    /// The scoring a bwa-mem3 command line aligns with, read from the words after
    /// its `mem` subcommand ([`mem_word_end`]; the whole command when it has
    /// none), up to the first `|`, `;`, `&&` or `||`:
    ///
    /// - `-A`, `-B`, `-O` and `-E` (value attached or as the next word; `-O`/`-E`
    ///   take `del[,ins]`). As in bwa, `-A` scales the defaults it does not
    ///   override, unless an `-x` preset is given (bwa then skips that scaling).
    /// - `-x intractg|pacbio|pbref|ont2d`, whose presets set `-B`/`-O`/`-E` (as
    ///   absolute values) where those options are not given.
    /// - `--meth` / `--meth=emseq|em-seq|bisulfite|taps` and `--meth-scoring
    ///   collapsed|genomic|neutral`: bisulfite alignment with collapsed scoring
    ///   (the default, except for `--meth=taps`, whose default is neutral) sets
    ///   `-B 2` (scaled by `-A`) unless `-B` is given.
    ///
    /// A value that does not parse is ignored, as is an unknown `-x` mode. Words
    /// are split on whitespace without interpreting shell quoting, so an
    /// option-like word inside a quoted argument (such as `-R '... -B 9 ...'`) is
    /// read as an option; `fgumi realign --help` documents this.
    #[must_use]
    pub(crate) fn from_command(command: &str) -> Self {
        let words: Vec<&str> = mem_option_words(command).collect();
        let mut set = ScoringWords::default();
        let mut at = 0;
        while at < words.len() {
            let word = words[at];
            if let Some(chem) = word.strip_prefix("--meth=") {
                set.meth_taps = Some(chem == "taps");
            } else if word == "--meth" {
                set.meth_taps = Some(false);
            } else if let Some(rest) = word.strip_prefix("--meth-scoring") {
                let mode = match rest.strip_prefix('=') {
                    Some(mode) => mode,
                    None if rest.is_empty() => {
                        at += 1;
                        words.get(at).copied().unwrap_or_default()
                    }
                    None => "",
                };
                set.meth_scoring = match mode {
                    "collapsed" => Some(MethScoring::Collapsed),
                    "genomic" => Some(MethScoring::Genomic),
                    "neutral" => Some(MethScoring::Neutral),
                    _ => set.meth_scoring,
                };
            } else if let Some(flag @ ('A' | 'B' | 'O' | 'E' | 'x')) =
                word.strip_prefix('-').and_then(|f| f.chars().next())
            {
                // The value is attached to the option, or else the next word.
                let attached = &word[2..];
                let value = if attached.is_empty() {
                    at += 1;
                    words.get(at).copied().unwrap_or_default()
                } else {
                    attached
                };
                match flag {
                    'A' => set.a = value.parse().ok().or(set.a),
                    'B' => set.b = value.parse().ok().or(set.b),
                    'O' => set.open = parse_pair(value).or(set.open),
                    'E' => set.extend = parse_pair(value).or(set.extend),
                    _ => set.mode = Some(value),
                }
            }
            at += 1;
        }
        Self::resolve(&set)
    }

    /// Applies bwa-mem3's defaults to what a command set, in bwa-mem3's order:
    /// the `-x` preset (or, without one, `-A` scaling), then `--meth`'s defaults.
    fn resolve(set: &ScoringWords<'_>) -> Self {
        let d = Self::default();
        let a = set.a.unwrap_or(d.a);
        // bwa-mem3 scales the defaults by `-A` only without an `-x` preset
        // (`if (mode) { preset } else update_a(..)`): under a preset, the
        // defaults the preset does not set stay unscaled.
        let preset = matches!(set.mode, Some("intractg" | "pacbio" | "pbref" | "ont2d"));
        let scale = if preset { 1 } else { a };
        let (mut o_del, mut o_ins) = set.open.unwrap_or((d.o_del * scale, d.o_ins * scale));
        let (mut e_del, mut e_ins) = set.extend.unwrap_or((d.e_del * scale, d.e_ins * scale));
        let mut b = set.b.unwrap_or(d.b * scale);
        match set.mode {
            Some("intractg") => {
                if set.open.is_none() {
                    (o_del, o_ins) = (16, 16);
                }
                if set.b.is_none() {
                    b = 9;
                }
            }
            Some("pacbio" | "pbref" | "ont2d") => {
                if set.open.is_none() {
                    (o_del, o_ins) = (1, 1);
                }
                if set.extend.is_none() {
                    (e_del, e_ins) = (1, 1);
                }
                if set.b.is_none() {
                    b = 1;
                }
            }
            _ => {}
        }
        if let Some(taps) = set.meth_taps {
            let default_mode = if taps { MethScoring::Neutral } else { MethScoring::Collapsed };
            if set.meth_scoring.unwrap_or(default_mode) == MethScoring::Collapsed && set.b.is_none()
            {
                b = Self::METH_MISMATCH * a;
            }
        }
        Self { a, b, o_del, e_del, o_ins, e_ins }
    }

    /// Scores an alignment from its BAM CIGAR ops and the number of mismatching
    /// aligned bases: `a` per matching aligned base, `-b` per mismatch, and
    /// `-(o + e * len)` per insertion or deletion run. Clipped bases score 0, and
    /// so does an alignment whose total is negative: bwa never reports a
    /// negative `AS`.
    #[must_use]
    pub(crate) fn score(&self, ops: &[u32], mismatches: u32) -> i32 {
        let mut aligned: i64 = 0;
        let mut penalty: i64 = 0;
        for &op in ops {
            let len = i64::from(op >> 4);
            match op & 0xF {
                0 | 7 | 8 => aligned += len,
                1 => penalty += i64::from(self.o_ins) + i64::from(self.e_ins) * len,
                2 => penalty += i64::from(self.o_del) + i64::from(self.e_del) * len,
                _ => {}
            }
        }
        let mismatches = i64::from(mismatches).min(aligned);
        let score =
            (aligned - mismatches) * i64::from(self.a) - mismatches * i64::from(self.b) - penalty;
        i32::try_from(score.max(0)).unwrap_or(i32::MAX)
    }
}

/// The words of `command` that bwa reads as its `mem` options: those after its
/// `mem` word ([`mem_word_end`]; the whole command when it has none), up to the
/// first `|`, `;`, `&&` or `||`. They are split on whitespace without
/// interpreting shell quoting, so an option-like word inside a quoted argument
/// counts as an option. The one rule every reader of a command's options uses.
pub(crate) fn mem_option_words(command: &str) -> impl Iterator<Item = &str> {
    let start = mem_word_end(command).unwrap_or(0);
    command[start..].split_whitespace().take_while(|w| !matches!(*w, "|" | ";" | "&&" | "||"))
}

/// Whether `command` passes bwa-mem3's `--meth` (`--meth` or `--meth=<chemistry>`)
/// among its [`mem_option_words`], as [`BwaScoring::from_command`] reads it.
pub(crate) fn command_passes_meth(command: &str) -> bool {
    mem_option_words(command).any(|w| w == "--meth" || w.starts_with("--meth="))
}

/// The byte offset just past the first `mem` word of `command` — `mem` with
/// whitespace (or the start of the command) before it and whitespace (or the end)
/// after it — or `None` when there is none. The one rule both reading a
/// command's options ([`BwaScoring::from_command`]) and inserting one
/// ([`super::ResolvedAligner::with_xa_max_hits`]) use to find bwa's subcommand.
#[must_use]
pub(crate) fn mem_word_end(command: &str) -> Option<usize> {
    let mut search = 0;
    while let Some(found) = command[search..].find("mem") {
        let begin = search + found;
        let end = begin + "mem".len();
        let before_ok = command[..begin].chars().next_back().is_none_or(char::is_whitespace);
        let after_ok = command[end..].chars().next().is_none_or(char::is_whitespace);
        if before_ok && after_ok {
            return Some(end);
        }
        search = end;
    }
    None
}

/// Parses bwa's `del[,ins]` gap value (`5` sets both).
fn parse_pair(value: &str) -> Option<(i32, i32)> {
    match value.split_once(',') {
        Some((del, ins)) => Some((del.parse().ok()?, ins.parse().ok()?)),
        None => value.parse().ok().map(|v| (v, v)),
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    const D: BwaScoring = BwaScoring { a: 1, b: 4, o_del: 6, e_del: 1, o_ins: 6, e_ins: 1 };

    /// The scoring read from a bwa-mem3 command line.
    #[rstest]
    #[case::none("bwa-mem3 mem -t 4 -p ref.fa /dev/stdin", D)]
    #[case::mismatch("bwa-mem3 mem -B 6 ref.fa -", BwaScoring { b: 6, ..D })]
    #[case::attached_values(
        "bwa mem -A2 -B3 -O5,7 -E2 ref.fa -",
        BwaScoring { a: 2, b: 3, o_del: 5, e_del: 2, o_ins: 7, e_ins: 2 }
    )]
    #[case::match_scales_defaults(
        "bwa mem -A 2 ref.fa -",
        BwaScoring { a: 2, b: 8, o_del: 12, e_del: 2, o_ins: 12, e_ins: 2 }
    )]
    #[case::single_gap_values(
        "bwa mem -O 5 -E 3 ref.fa -",
        BwaScoring { o_del: 5, e_del: 3, o_ins: 5, e_ins: 3, ..D }
    )]
    #[case::stops_at_pipe("bwa mem -B 6 ref.fa - | other -B 9", BwaScoring { b: 6, ..D })]
    #[case::only_after_mem("env -B 9 bwa mem -B 5 ref.fa -", BwaScoring { b: 5, ..D })]
    #[case::tab_before_mem("bwa\tmem -B 6 ref.fa -", BwaScoring { b: 6, ..D })]
    #[case::no_mem_word("my-aligner -B 5 {ref}", BwaScoring { b: 5, ..D })]
    #[case::unparsable_ignored("bwa mem -B x ref.fa -", D)]
    #[case::meth_defaults("bwa-mem3 mem --meth ref.fa -", BwaScoring { b: 2, ..D })]
    #[case::meth_emseq_chemistry("bwa-mem3 mem --meth=emseq ref.fa -", BwaScoring { b: 2, ..D })]
    #[case::meth_scaled_by_match(
        "bwa-mem3 mem --meth -A 2 ref.fa -",
        BwaScoring { a: 2, b: 4, o_del: 12, e_del: 2, o_ins: 12, e_ins: 2 }
    )]
    #[case::meth_explicit_mismatch("bwa-mem3 mem --meth -B 3 ref.fa -", BwaScoring { b: 3, ..D })]
    #[case::meth_genomic("bwa-mem3 mem --meth --meth-scoring genomic ref.fa -", D)]
    #[case::meth_neutral_attached("bwa-mem3 mem --meth --meth-scoring=neutral ref.fa -", D)]
    #[case::meth_collapsed_explicit(
        "bwa-mem3 mem --meth --meth-scoring collapsed ref.fa -",
        BwaScoring { b: 2, ..D }
    )]
    #[case::taps_defaults_to_neutral("bwa-mem3 mem --meth=taps ref.fa -", D)]
    #[case::taps_collapsed(
        "bwa-mem3 mem --meth=taps --meth-scoring collapsed ref.fa -",
        BwaScoring { b: 2, ..D }
    )]
    #[case::meth_scoring_without_meth("bwa-mem3 mem --meth-scoring collapsed ref.fa -", D)]
    #[case::preset_intractg(
        "bwa-mem3 mem -x intractg ref.fa -",
        BwaScoring { b: 9, o_del: 16, o_ins: 16, ..D }
    )]
    #[case::preset_pacbio_attached(
        "bwa-mem3 mem -xpacbio ref.fa -",
        BwaScoring { a: 1, b: 1, o_del: 1, e_del: 1, o_ins: 1, e_ins: 1 }
    )]
    #[case::preset_ont2d_explicit_wins(
        "bwa-mem3 mem -x ont2d -B 3 -O 2 ref.fa -",
        BwaScoring { a: 1, b: 3, o_del: 2, e_del: 1, o_ins: 2, e_ins: 1 }
    )]
    #[case::preset_unknown("bwa-mem3 mem -x short ref.fa -", D)]
    // Under a preset bwa does not scale by `-A`: intractg leaves `-E` at 1.
    #[case::preset_intractg_with_match_score(
        "bwa-mem3 mem -x intractg -A 2 ref.fa -",
        BwaScoring { a: 2, b: 9, o_del: 16, e_del: 1, o_ins: 16, e_ins: 1 }
    )]
    fn from_command_cases(#[case] command: &str, #[case] expected: BwaScoring) {
        assert_eq!(BwaScoring::from_command(command), expected);
    }

    /// Scored from the CIGAR and the mismatch count; never negative, as bwa never
    /// reports a negative `AS`.
    #[rstest]
    #[case::all_match("20M", 0, 20)]
    #[case::one_mismatch("20M", 1, 15)]
    #[case::soft_clip("2S18M", 0, 18)]
    #[case::insertion("10M1I9M", 0, 12)]
    #[case::deletion("10M2D10M", 0, 12)]
    #[case::mismatch_and_deletion("10M2D10M", 2, 2)]
    #[case::negative_clamps_to_zero("20M", 16, 0)]
    fn score_cases(#[case] cigar: &str, #[case] mismatches: u32, #[case] expected: i32) {
        let ops = fgumi_raw_bam::parse_cigar_string_ops(cigar).expect("cigar");
        assert_eq!(D.score(&ops, mismatches), expected);
    }

    #[rstest]
    #[case::subcommand("bwa-mem3 mem -p ref.fa -", Some(12))]
    #[case::first_word("mem -p ref.fa -", Some(3))]
    #[case::tab("bwa\tmem -p", Some(7))]
    #[case::at_end("bwa mem", Some(7))]
    #[case::path_component("/opt/mem/bwa mem x", Some(16))]
    #[case::prefix_word("tool memo x", None)]
    #[case::suffix_word("tool xmem x", None)]
    #[case::none("minimap2 -ax sr ref.fa -", None)]
    fn mem_word_end_cases(#[case] command: &str, #[case] expected: Option<usize>) {
        assert_eq!(mem_word_end(command), expected);
    }

    /// `--meth` is read from the same words as the scoring options: after `mem`,
    /// and only up to the first shell separator.
    #[rstest]
    #[case::flag("bwa-mem3 mem --meth -p ref.fa -", true)]
    #[case::chemistry("bwa-mem3 mem --meth=taps -p ref.fa -", true)]
    #[case::absent("bwa-mem3 mem -p ref.fa -", false)]
    #[case::before_mem("env --meth bwa-mem3 mem -p ref.fa -", false)]
    #[case::after_pipe("bwa-mem3 mem -p ref.fa - | tool --meth", false)]
    #[case::after_and("bwa-mem3 mem -p ref.fa - && tool --meth", false)]
    fn command_passes_meth_cases(#[case] command: &str, #[case] expected: bool) {
        assert_eq!(command_passes_meth(command), expected);
    }
}
