use std::env;

fn main() {
    let semel_lib = env::var("SEMEL_LIB_PATH")
        .unwrap_or_else(|_| "/home/marcello/dev/semel/result/lib".to_string());

    println!("cargo:rustc-link-search=native={}", semel_lib);
    println!("cargo:rustc-link-lib=dylib=semel");

    // Propagate the library path so downstream crates can set rpath
    println!("cargo:semel_lib_path={}", semel_lib);

    // Re-run if the env var changes
    println!("cargo:rerun-if-env-changed=SEMEL_LIB_PATH");
}
