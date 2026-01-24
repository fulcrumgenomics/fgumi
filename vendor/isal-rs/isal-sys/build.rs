use std::path::PathBuf;
use std::{
    io::{self, Write},
    process::{Command, Stdio},
};

fn main() {
    println!("cargo:rerun-if-changed=build.rs");

    let is_static = cfg!(feature = "static");
    let is_shared = cfg!(feature = "shared");
    let rust_target = std::env::var("TARGET").unwrap();
    let target = if rust_target.starts_with("riscv64") {
        let (_cpu, rest) = rust_target.split_once('-').unwrap();
        format!("riscv64-{rest}")
    } else {
        rust_target
    };
    let out_dir = PathBuf::from(&std::env::var("OUT_DIR").unwrap());

    // Copy isa-l source into out; not allow to modify things outside of out dir
    let src_dir = out_dir.join("isa-l");
    if src_dir.exists() {
        std::fs::remove_dir_all(&src_dir).unwrap(); // maybe from a previous build
    }
    copy_dir::copy_dir("isa-l", &src_dir).unwrap();

    let install_path = std::env::var("ISAL_INSTALL_PREFIX")
        .map(|p| PathBuf::from(&p).clone())
        .unwrap_or(out_dir.clone());

    let current_dir = std::env::current_dir().unwrap();
    std::env::set_current_dir(&src_dir).unwrap();

    // build from source
    #[cfg(not(feature = "use-system-isal"))]
    {
        #[cfg(not(target_os = "windows"))]
        let cmd = {
            let status = Command::new("./autogen.sh")
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .output()
                .unwrap();
            io::stdout().write_all(&status.stdout).unwrap();
            io::stderr().write_all(&status.stderr).unwrap();
            if !status.status.success() {
                panic!("autogen failed");
            }

            let compiler = cc::Build::new().get_compiler();
            let cflags = compiler.cflags_env().into_string().unwrap();

            let mut configure_args = vec![
                format!("--prefix={}", install_path.display()),
                format!("--host={}", target),
                format!("--enable-static={}", if is_static { "yes" } else { "no" }),
                format!("--enable-shared={}", if is_shared { "yes" } else { "no" }),
                format!("CFLAGS={}", cflags),
                format!("CC={}", compiler.path().display()),
            ];

            if !cfg!(target_os = "macos") {
                let ldflag = if is_static { "static" } else { "shared" };
                configure_args.push(format!("LDFLAGS=-{}", ldflag));
                configure_args.push("--with-pic=yes".to_string());
            }

            let status = Command::new("./configure")
                .args(&configure_args)
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .output()
                .unwrap();
            io::stdout().write_all(&status.stdout).unwrap();
            io::stderr().write_all(&status.stderr).unwrap();
            if !status.status.success() {
                panic!("configure failed");
            }

            Command::new("make")
                .args(&["install-libLTLIBRARIES"])
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
        };

        #[cfg(target_os = "windows")]
        let mut cmd = {
            Command::new("nmake")
                .args(["-f", "Makefile.nmake"])
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
        };

        std::env::set_current_dir(&current_dir).unwrap();

        let output = cmd.unwrap().wait_with_output().unwrap();
        io::stdout().write_all(&output.stdout).unwrap();
        io::stderr().write_all(&output.stderr).unwrap();
        if !output.status.success() {
            panic!("Building isa-l failed");
        }
    }

    let libname = if cfg!(target_os = "windows") {
        println!("cargo:rustc-link-search=native={}", src_dir.display());
        if cfg!(feature = "static") {
            "isa-l_static"
        } else {
            "isa-l"
        }
    } else {
        for subdir in ["bin", "lib", "lib64"] {
            let search_path = install_path.join(subdir);
            println!("cargo:rustc-link-search=native={}", search_path.display());
        }
        "isal"
    };

    #[cfg(feature = "static")]
    println!("cargo:rustc-link-lib=static={}", libname);

    #[cfg(feature = "shared")]
    println!("cargo:rustc-link-lib={}", libname);

    #[cfg(feature = "regenerate-bindings")]
    {
        let out = PathBuf::from(&(format!("{}/igzip_lib.rs", std::env::var("OUT_DIR").unwrap())));
        bindgen::Builder::default()
            // The input header we would like to generate
            // bindings for.
            .header("isa-l/include/igzip_lib.h")
            // Tell cargo to invalidate the built crate whenever any of the
            // included header files changed.
            .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
            .blocklist_type("__uint64_t_")
            .blocklist_type("__size_t")
            // Replaced by libc::FILE
            .blocklist_type("FILE")
            .blocklist_type("_IO_FILE")
            .blocklist_type("_IO_codecvt")
            .blocklist_type("_IO_wide_data")
            .blocklist_type("_IO_marker")
            .blocklist_type("_IO_lock_t")
            // Replaced by i64
            .blocklist_type("LARGE_INTEGER")
            // Replaced by libc::timespec
            .blocklist_type("timespec")
            // etc
            .blocklist_type("__time_t")
            .blocklist_type("__syscall_slong_t")
            // .blocklist_type("__off64_t")
            .blocklist_type("__off_t")
            .size_t_is_usize(true)
            .generate()
            .expect("Unable to generate bindings")
            .write_to_file(out)
            .unwrap();
    }
}
