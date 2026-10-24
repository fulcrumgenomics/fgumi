//! R guard shared by the duplex- and simplex-metrics PDF tests.

use std::process::Command;

/// Whether the `Rscript` on `PATH` can leave the creation and modification dates out of a PDF,
/// through the `timestamp` argument R 4.5.0 added to `pdf()`.
pub fn r_pdf_omits_dates() -> bool {
    Command::new("Rscript")
        .args(["-e", r#"stopifnot("timestamp" %in% names(formals(grDevices::pdf)))"#])
        .output()
        .is_ok_and(|output| output.status.success())
}
