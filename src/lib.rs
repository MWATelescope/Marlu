// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Core code to describe coordinate transformations, Jones matrices, etc.

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
pub mod sexagesimal;

#[cfg(feature = "io")]
pub mod io;

#[cfg(feature = "cuda")]
pub mod cuda;

// Re-exports.
pub use jones::Jones;
pub use pos::{
    azel::AzEl,
    earth::{Ellipsoid, LatLngHeight},
    enh::ENH,
    hadec::HADec,
    lmn::LMN,
    pal, precession,
    radec::RADec,
    uvw::UVW,
    xyz::{XyzGeocentric, XyzGeodetic},
};

pub use context::MarluVisContext;
pub use erfa_sys;
pub use hifitime;
pub use ndarray;
pub use num_complex;
pub use num_complex::Complex;
pub use num_traits;
pub use rayon;

// If "mwalib" is enabled, re-export the crate here, as well its re-exported
// crates.
cfg_if::cfg_if! {
    if #[cfg(feature = "mwalib")] {
        pub use mwalib;
        pub use mwalib::{fitsio, fitsio_sys};
    }
}

// If "io" is enabled, re-export rubbl_casatables here.
#[cfg(feature = "io")]
pub use rubbl_casatables;

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
