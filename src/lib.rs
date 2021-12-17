// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Core code to describe coordinate transformations, Jones matrices, etc.

#[allow(non_camel_case_types)]
pub type c32 = num_complex::Complex<f32>;
#[allow(non_camel_case_types)]
pub type c64 = num_complex::Complex<f64>;

pub mod constants;
pub mod jones;
pub mod math;
pub mod pos;
pub mod sexagesimal;
pub mod time;

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
