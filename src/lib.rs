// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Core code to describe coordinate transformations, Jones matrices, etc.

#![deny(clippy::all)]
#![warn(clippy::missing_safety_doc)]
// #![warn(clippy::missing_errors_doc)]
#![warn(clippy::float_cmp)]
#![warn(clippy::if_not_else)]
#![warn(clippy::explicit_iter_loop)]
#![warn(clippy::wildcard_imports)]
#![warn(clippy::cloned_instead_of_copied)]
#![warn(clippy::cognitive_complexity)]
#![warn(clippy::type_repetition_in_bounds)]
#![warn(clippy::redundant_closure_for_method_calls)]
#![warn(clippy::doc_markdown)]
#![warn(clippy::match_same_arms)]
#![warn(clippy::default_trait_access)]
#![warn(clippy::semicolon_if_nothing_returned)]
#![warn(clippy::explicit_into_iter_loop)]
#![warn(clippy::inefficient_to_string)]
#![warn(clippy::needless_pass_by_value)]
#![warn(clippy::used_underscore_binding)]

#[allow(non_camel_case_types)]
pub type c32 = num_complex::Complex<f32>;
#[allow(non_camel_case_types)]
pub type c64 = num_complex::Complex<f64>;

pub mod averaging;
pub mod constants;
pub mod context;
pub mod jones;
pub mod math;
pub mod pos;
pub mod selection;
pub mod sexagesimal;

#[cfg(feature = "io")]
pub mod io;

#[cfg(feature = "cuda")]
pub mod cuda;

// Re-exports.
pub use context::{History, MwaObsContext, ObsContext, VisContext};
pub use jones::Jones;
pub use pos::{
    azel::AzEl,
    earth::{Ellipsoid, LatLngHeight},
    enh::ENH,
    hadec::HADec,
    lmn::{LmnRime, LMN},
    pal, precession,
    radec::RADec,
    uvw::UVW,
    xyz::{XyzGeocentric, XyzGeodetic},
};
pub use selection::{SelectionError, VisSelection};

pub use erfa_sys;
pub use hifitime;
pub use ndarray;
pub use num_complex;
pub use num_complex::Complex;
pub use num_traits;
pub use rayon;

// Include the generated built.rs code into our library
pub mod built_info {
    // The file has been placed there by the build script.
    include!(concat!(env!("OUT_DIR"), "/built.rs"));
}

// If "mwalib" is enabled, re-export the crate here, as well its re-exported
// crates.
cfg_if::cfg_if! {
    if #[cfg(feature = "mwalib")] {
        pub use mwalib;
        pub use mwalib::{fitsio, fitsio_sys};
    }
}

// If "io" is enabled, re-export rubbl_casatables here.
cfg_if::cfg_if! {
    if #[cfg(feature = "io")] {
        pub use rubbl_casatables;
        pub use io::{MeasurementSetWriter, UvfitsWriter, VisWritable, UvfitsWriteError};
    }
}

// If "cuda" is enabled, re-export cuda-runtime-sys here.
#[cfg(feature = "cuda")]
pub use cuda_runtime_sys;

#[cfg(test)]
#[test]
fn hifitime_works_as_expected() {
    use hifitime::Epoch;

    let gps = 1065880128.0;
    let epoch = Epoch::from_gpst_seconds(gps);
    approx::assert_abs_diff_eq!(epoch.as_gpst_seconds(), gps);

    let jd_utc = 2444244.5;
    let epoch = Epoch::from_jde_utc(jd_utc);
    approx::assert_abs_diff_eq!(epoch.as_jde_utc_days(), jd_utc);
    approx::assert_abs_diff_eq!(epoch.as_gpst_seconds(), 0.0);
}
