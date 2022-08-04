// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// This code is adapted from pkg-config-rs
// (https://github.com/rust-lang/pkg-config-rs).
#[cfg(feature = "cuda")]
#[allow(clippy::if_same_then_else, clippy::needless_bool)]
fn infer_static(name: &str) -> bool {
    if std::env::var(&format!("{}_STATIC", name.to_uppercase())).is_ok() {
        true
    } else if std::env::var(&format!("{}_DYNAMIC", name.to_uppercase())).is_ok() {
        false
    } else if std::env::var("PKG_CONFIG_ALL_STATIC").is_ok() {
        true
    } else if std::env::var("PKG_CONFIG_ALL_DYNAMIC").is_ok() {
        false
    } else {
        false
    }
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");

    // Gather build time info
    built::write_built_file().expect("Failed to acquire build-time information");

    #[cfg(feature = "cuda")]
    {
        // Link CUDA. If the library path manually specified, search there.
        if let Ok(lib_dir) = std::env::var("CUDA_LIB") {
            println!("cargo:rustc-link-search=native={}", lib_dir);
        }

        if infer_static("cuda") {
            // CUDA ships its static library as cudart_static.a, not cudart.a
            println!("cargo:rustc-link-lib=static=cudart_static");
        } else {
            println!("cargo:rustc-link-lib=cudart");
        }

        #[cfg(feature = "cuda-static")]
        println!("cargo:rustc-link-lib=static=cudart_static");
    }
}
