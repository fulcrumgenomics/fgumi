//! Shared serde helpers for writing `f64` metric fields in fgbio's `Metric` TSV format.
//!
//! fgbio serializes doubles with `Double.toString`, which renders non-finite values
//! as `Infinity`, `-Infinity`, and `NaN` — the exact tokens `Double.parseDouble`
//! accepts on read. Rust's default `f64` serialization (via ryu) instead emits
//! `inf`/`-inf`/`NaN`, and `inf`/`-inf` are **unparseable** by fgbio's `Metric.read`.
//!
//! Every `f64` field in an fgumi metric that is meant to be read back by fgbio must
//! carry `#[serde(with = "crate::float")]` so non-finite values round-trip through
//! both tools. Finite whole numbers are written without a trailing `.0` (e.g. `0`
//! rather than `0.0`), matching fgbio's `DecimalFormat` output for integral values;
//! other finite values keep full precision (numeric equality — not byte-for-byte
//! formatting — is the metric contract, see the crate `format_float` docs).

use serde::{Deserialize, Deserializer, Serializer};

/// Serializes an `f64` in fgbio's `Metric` TSV format.
///
/// Non-finite values become `Infinity`/`-Infinity`/`NaN` (readable by fgbio's
/// `Double.parseDouble`); finite whole numbers drop the fractional part; all other
/// finite values keep full precision.
///
/// Intended for use as `#[serde(with = "fgumi_metrics::float")]` on `f64` metric
/// fields so every crate's metric TSVs format floats consistently.
///
/// # Errors
///
/// Propagates any error from the underlying `serializer`.
#[allow(clippy::trivially_copy_pass_by_ref, reason = "serde requires &T signature")]
pub fn serialize<S>(value: &f64, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let value = *value;
    if value.is_nan() {
        return serializer.serialize_str("NaN");
    }
    if value.is_infinite() {
        return serializer.serialize_str(if value > 0.0 { "Infinity" } else { "-Infinity" });
    }

    // Render an integral value without a decimal point (fgbio's DecimalFormat does the same).
    #[expect(
        clippy::cast_precision_loss,
        reason = "i64::MAX boundary check, exact precision not needed"
    )]
    let max_safe = i64::MAX as f64;
    if value.fract() == 0.0 && value.abs() < max_safe {
        #[expect(clippy::cast_possible_truncation, reason = "guarded by abs check above")]
        let int_value = value as i64;
        serializer.serialize_str(&format!("{int_value}"))
    } else {
        serializer.serialize_str(&format!("{value}"))
    }
}

/// Deserializes an `f64` written by [`serialize`], accepting the non-finite
/// tokens `NaN`, `Infinity`, and `-Infinity` in addition to ordinary numbers.
///
/// # Errors
///
/// Returns an error if the underlying `deserializer` fails or the value is
/// neither a recognized non-finite token nor a parseable number.
pub fn deserialize<'de, D>(deserializer: D) -> Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    let s: String = String::deserialize(deserializer)?;
    match s.as_str() {
        "NaN" => Ok(f64::NAN),
        "Infinity" => Ok(f64::INFINITY),
        "-Infinity" => Ok(f64::NEG_INFINITY),
        _ => s.parse().map_err(serde::de::Error::custom),
    }
}

#[cfg(test)]
mod tests {
    use fgoxide::io::DelimFile;
    use proptest::prelude::*;
    use rstest::rstest;
    use serde::{Deserialize, Serialize};
    use tempfile::NamedTempFile;

    #[derive(Debug, Serialize, Deserialize)]
    struct Row {
        #[serde(with = "super")]
        value: f64,
    }

    /// Writes a one-column TSV via the real `DelimFile` path and returns the value cell
    /// (the second physical line, after the `value` header).
    fn to_cell(value: f64) -> String {
        let file = NamedTempFile::new().expect("temp file");
        DelimFile::default().write_tsv(file.path(), [Row { value }]).expect("write");
        let content = std::fs::read_to_string(file.path()).expect("read");
        content.lines().nth(1).expect("value row").to_string()
    }

    /// Writes `value` and reads the single row back through the real `DelimFile` path,
    /// returning the round-tripped `f64`.
    fn round_trip(value: f64) -> f64 {
        let file = NamedTempFile::new().expect("temp file");
        DelimFile::default().write_tsv(file.path(), [Row { value }]).expect("write");
        let rows: Vec<Row> = DelimFile::default().read_tsv(file.path()).expect("read");
        assert_eq!(rows.len(), 1);
        rows[0].value
    }

    /// Each value serializes to fgbio's exact TSV token: non-finite values become
    /// `Infinity`/`-Infinity`/`NaN`, integral values drop the decimal point, and
    /// fractional values keep full precision.
    #[rstest]
    #[case::pos_infinity(f64::INFINITY, "Infinity")]
    #[case::neg_infinity(f64::NEG_INFINITY, "-Infinity")]
    #[case::nan(f64::NAN, "NaN")]
    #[case::zero(0.0, "0")]
    #[case::positive_integral(42.0, "42")]
    #[case::negative_integral(-7.0, "-7")]
    #[case::half(0.5, "0.5")]
    #[case::fractional(1.25, "1.25")]
    fn writes_fgbio_token(#[case] value: f64, #[case] expected: &str) {
        assert_eq!(to_cell(value), expected);
    }

    /// Each signed infinity round-trips back to an infinity of the same sign.
    #[rstest]
    #[case::pos_infinity(f64::INFINITY)]
    #[case::neg_infinity(f64::NEG_INFINITY)]
    fn round_trips_signed_infinity(#[case] value: f64) {
        let read = round_trip(value);
        assert!(read.is_infinite() && read.signum() == value.signum());
    }

    #[test]
    fn round_trips_nan() {
        assert!(round_trip(f64::NAN).is_nan());
    }

    proptest! {
        /// Any finite value round-trips through the write/read path to an equal `f64`
        /// (serialization uses Rust's shortest round-trippable representation).
        #[test]
        fn round_trips_finite(value in -1e15_f64..1e15_f64) {
            let read = round_trip(value);
            prop_assert!((read - value).abs() < f64::EPSILON, "round-trip {value}");
        }
    }
}
