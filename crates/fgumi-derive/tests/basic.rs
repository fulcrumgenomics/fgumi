use fgumi_derive::multi_options;

#[multi_options("test", "Test Options")]
#[derive(clap::Args, Debug, Clone)]
pub struct TestOptions {
    /// A required string field.
    #[arg(long)]
    pub name: String,

    /// An optional path.
    #[arg(long)]
    pub output: Option<std::path::PathBuf>,

    /// A field with a default.
    #[arg(long, default_value_t = 42)]
    pub count: usize,
}

impl Default for TestOptions {
    fn default() -> Self {
        Self { name: String::new(), output: None, count: 42 }
    }
}

#[test]
fn test_multi_options_generates_struct() {
    let multi = MultiTestOptions {
        test_name: Some("hello".to_string()),
        test_output: None,
        test_count: 42,
    };
    let original = multi.validate().unwrap();
    assert_eq!(original.name, "hello");
    assert_eq!(original.output, None);
    assert_eq!(original.count, 42);
}

#[test]
fn test_multi_options_validate_rejects_missing_required() {
    let multi = MultiTestOptions { test_name: None, test_output: None, test_count: 42 };
    assert!(multi.validate().is_err());
}

#[test]
fn test_multi_options_non_default_value_preserved() {
    let multi = MultiTestOptions {
        test_name: Some("world".to_string()),
        test_output: Some(std::path::PathBuf::from("/tmp/out")),
        test_count: 7,
    };
    let original = multi.validate().unwrap();
    assert_eq!(original.name, "world");
    assert_eq!(original.output, Some(std::path::PathBuf::from("/tmp/out")));
    assert_eq!(original.count, 7);
}

#[test]
fn test_multi_options_from_original() {
    let original = TestOptions { name: "test".to_string(), output: None, count: 99 };
    let multi = MultiTestOptions::from(original);
    assert_eq!(multi.test_name, Some("test".to_string()));
    assert_eq!(multi.test_output, None);
    assert_eq!(multi.test_count, 99);
}

#[test]
fn test_multi_options_roundtrip() {
    let original = TestOptions {
        name: "roundtrip".to_string(),
        output: Some(std::path::PathBuf::from("/x")),
        count: 123,
    };
    let multi = MultiTestOptions::from(original);
    let back = multi.validate().unwrap();
    assert_eq!(back.name, "roundtrip");
    assert_eq!(back.output, Some(std::path::PathBuf::from("/x")));
    assert_eq!(back.count, 123);
}
