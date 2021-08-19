// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Useful constants.

use std::f64::consts::PI;

/// Speed of light \[metres/second\]
pub const VEL_C: f64 = erfa_sys::ERFA_CMPS;

/// Seconds per day (86400)
pub const DAYSEC: f64 = erfa_sys::ERFA_DAYSEC;
/// Seconds of time to radians (7.272205216643039903848712e-5).
pub const DS2R: f64 = erfa_sys::ERFA_DS2R;
/// Hour angle to radians (15 / 180 * PI).
pub const DH2R: f64 = 15.0 / 180.0 * PI;
/// Ratio of a solar day to a sidereal day (24/23.9344696 = 1.002737909).
pub const SOLAR2SIDEREAL: f64 = 24.0 / 23.9344696;

/// MWA latitude \[radians\]
pub const MWA_LAT_RAD: f64 = -0.4660608448386394;
/// MWA latitude \[degrees\]
pub const MWA_LAT_DEG: f64 = MWA_LAT_RAD * 180.0 / PI;
/// MWA longitude \[radians\]
pub const MWA_LONG_RAD: f64 = 2.0362898668561042;
/// MWA latitude \[degrees\]
pub const MWA_LONG_DEG: f64 = MWA_LONG_RAD * 180.0 / PI;

/// MWA height (a.k.a. altitude) \[metres\]
pub const MWA_HEIGHT_M: f64 = 377.827;

// cotter's constants. Useful for being more precise when converting geocentric
// XYZ to geodetic XYZ!
/// cotter's MWA latitude on Earth in radians. Use [MWA_LAT_RAD] unless you know
/// what you're doing.
pub const COTTER_MWA_LATITUDE_RADIANS: f64 = -0.46606083776035967;
/// cotter's MWA longitude on Earth in radians. Use [MWA_LONG_RAD] unless you
/// know what you're doing.
pub const COTTER_MWA_LONGITUDE_RADIANS: f64 = 2.0362897754687257;
/// cotter's MWA altitude in metres. Use [MWA_HEIGHT_M] unless you know what
/// you're doing.
pub const COTTER_MWA_HEIGHT_METRES: f64 = 377.0;

/// This is the number of seconds from 1900 Jan 1 and 1980 Jan 5. The GPS epoch
/// is 1980 Jan 5, but `hifitime` uses 1900 for everything; subtracting this
/// number from the result of `hifitime::Epoch::as_gpst_seconds` gives the
/// expected GPS time.
pub const HIFITIME_GPS_FACTOR: f64 =
    hifitime::SECONDS_PER_YEAR * 80.0 + hifitime::SECONDS_PER_DAY * 4.0;

/// The number of seconds between 1858-11-17T00:00:00 (MJD epoch, used by
/// casacore) and 1900-01-01T00:00:00 (TAI epoch) is 1297728000. I'm using the
/// TAI epoch because that's well supported by `hifitime`, and `hifitime`
/// converts an epoch to many formats including JD, and accounts for leap
/// seconds.
pub const MJD_TAI_EPOCH_DIFF: f64 = 1297728000.0;

// Double check that our constants match mwalib's. This crate would just use
// mwalib's constants, but mwalib is an optional dependency.
#[cfg(feature = "mwalib")]
extern crate static_assertions as sa;
#[cfg(feature = "mwalib")]
sa::const_assert!(MWA_LAT_RAD == mwalib::MWA_LATITUDE_RADIANS);
#[cfg(feature = "mwalib")]
sa::const_assert!(MWA_LONG_RAD == mwalib::MWA_LONGITUDE_RADIANS);
#[cfg(feature = "mwalib")]
sa::const_assert!(MWA_HEIGHT_M == mwalib::MWA_ALTITUDE_METRES);
