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
#[expect(clippy::trivially_copy_pass_by_ref, reason = "serde requires &T signature")]
pub(crate) fn serialize<S>(value: &f64, serializer: S) -> Result<S::Ok, S::Error>
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
pub(crate) fn deserialize<'de, D>(deserializer: D) -> Result<f64, D::Error>
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

    #[test]
    fn writes_fgbio_readable_non_finite_tokens() {
        assert_eq!(to_cell(f64::INFINITY), "Infinity");
        assert_eq!(to_cell(f64::NEG_INFINITY), "-Infinity");
        assert_eq!(to_cell(f64::NAN), "NaN");
    }

    #[test]
    fn writes_integral_values_without_decimal_point() {
        assert_eq!(to_cell(0.0), "0");
        assert_eq!(to_cell(42.0), "42");
        assert_eq!(to_cell(-7.0), "-7");
    }

    #[test]
    fn keeps_full_precision_for_fractional_values() {
        assert_eq!(to_cell(0.5), "0.5");
        assert_eq!(to_cell(1.25), "1.25");
    }

    #[test]
    fn round_trips_non_finite_and_finite() {
        for &v in &[0.0_f64, 1.5, -3.25, 1000.0, f64::INFINITY, f64::NEG_INFINITY] {
            let file = NamedTempFile::new().expect("temp file");
            DelimFile::default().write_tsv(file.path(), [Row { value: v }]).expect("write");
            let rows: Vec<Row> = DelimFile::default().read_tsv(file.path()).expect("read");
            assert_eq!(rows.len(), 1);
            if v.is_infinite() {
                assert!(rows[0].value.is_infinite() && rows[0].value.signum() == v.signum());
            } else {
                assert!((rows[0].value - v).abs() < f64::EPSILON, "round-trip {v}");
            }
        }
    }
}
