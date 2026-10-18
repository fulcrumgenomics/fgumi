//! Helpers shared by this crate's unit tests.
//!
//! `cfg(test)` only, not `test-utils`: a dependent's tests may install their
//! own logger, and nothing here should be able to collide with it.

use std::sync::Once;

use parking_lot::Mutex;

/// One log record seen by the capture logger.
#[derive(Clone, Debug)]
pub(crate) struct CapturedRecord {
    /// The record's target (by default its module path).
    pub target: String,
    /// The record's level.
    pub level: log::Level,
    /// The formatted message.
    pub message: String,
}

/// Every record logged since the capture logger was installed.
static CAPTURED: Mutex<Vec<CapturedRecord>> = Mutex::new(Vec::new());

/// The capture logger: keeps every record, of every target and level.
struct Capture;

impl log::Log for Capture {
    fn enabled(&self, _metadata: &log::Metadata) -> bool {
        true
    }

    fn log(&self, record: &log::Record) {
        CAPTURED.lock().push(CapturedRecord {
            target: record.target().to_owned(),
            level: record.level(),
            message: record.args().to_string(),
        });
    }

    fn flush(&self) {}
}

/// Install the process-wide capture logger, if it is not installed yet.
///
/// `log` admits one logger per process, so a test in this crate that needs to
/// see log records calls this rather than installing its own; a second
/// `log::set_logger` would fail, and the next test to try would panic. Call it
/// before the code under test logs: records logged earlier are not kept.
///
/// # Panics
///
/// Panics if a logger other than this one was already installed.
pub(crate) fn install_log_capture() {
    static LOGGER: Capture = Capture;
    static INSTALL: Once = Once::new();
    INSTALL.call_once(|| {
        log::set_logger(&LOGGER).expect(
            "another logger is installed; capture logs through `test_support::install_log_capture`",
        );
        log::set_max_level(log::LevelFilter::Trace);
    });
}

/// The captured records from `target` whose message contains `needle`.
///
/// Tests share one buffer (`cargo test` runs them concurrently in one
/// process), so pick a `needle` only the calling test logs, such as its own
/// step names.
pub(crate) fn captured_logs(target: &str, needle: &str) -> Vec<CapturedRecord> {
    install_log_capture();
    CAPTURED
        .lock()
        .iter()
        .filter(|r| r.target == target && r.message.contains(needle))
        .cloned()
        .collect()
}
