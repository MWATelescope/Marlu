// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Handle (hour angle, declination) coordinates.

use crate::{constants::MWA_LAT_RAD, AzEl, RADec};

/// A struct containing an Hour Angle and Declination. All units are in radians.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[allow(clippy::upper_case_acronyms)]
pub struct HADec {
    /// Hour angle \[radians\]
    pub ha: f64,
    /// Declination \[radians\]
    pub dec: f64,
}

impl HADec {
    /// Make a new [`HADec`] struct from values in radians.
    pub fn from_radians(ha: f64, dec: f64) -> HADec {
        Self { ha, dec }
    }

    /// Make a new [`HADec`] struct from values in degrees.
    pub fn from_degrees(ha: f64, dec: f64) -> HADec {
        Self {
            ha: ha.to_radians(),
            dec: dec.to_radians(),
        }
    }

    /// Make a new [`HADec`] struct from values in radians.
    #[deprecated = "use `HADec::from_radians` instead"]
    pub fn new(ha_rad: f64, dec_rad: f64) -> HADec {
        Self::from_radians(ha_rad, dec_rad)
    }

    /// Make a new [`HADec`] struct from values in degrees.
    #[deprecated = "use `HADec::from_degrees` instead"]
    pub fn new_degrees(ha_deg: f64, dec_deg: f64) -> HADec {
        Self::from_degrees(ha_deg, dec_deg)
    }

    /// Given a local sidereal time, make a new [`RADec`] struct from a [`HADec`].
    pub fn to_radec(self, lst_rad: f64) -> RADec {
        RADec {
            ra: lst_rad - self.ha,
            dec: self.dec,
        }
    }

    /// Given a local sidereal time, make a new [`HADec`] struct from a [`RADec`].
    pub fn from_radec(radec: RADec, lst_rad: f64) -> HADec {
        Self {
            ha: lst_rad - radec.ra,
            dec: radec.dec,
        }
    }

    /// Convert the equatorial coordinates to horizon coordinates (azimuth and
    /// elevation), given the local latitude on Earth.
    ///
    /// Uses ERFA.
    pub fn to_azel(self, latitude_rad: f64) -> AzEl {
        let mut az = 0.0;
        let mut el = 0.0;
        unsafe { erfa_sys::eraHd2ae(self.ha, self.dec, latitude_rad, &mut az, &mut el) }
        AzEl::from_radians(az, el)
    }

    /// Convert the equatorial coordinates to horizon coordinates (azimuth and
    /// elevation) for the MWA's location.
    ///
    /// Uses ERFA.
    pub fn to_azel_mwa(self) -> AzEl {
        self.to_azel(MWA_LAT_RAD)
    }

    /// Calculate the distance between two sets of coordinates.
    ///
    /// Uses ERFA.
    pub fn separation(self, b: Self) -> f64 {
        unsafe { erfa_sys::eraSeps(self.ha, self.dec, b.ha, b.dec) }
    }

    /// Get the [parallactic
    /// angle](https://en.wikipedia.org/wiki/Parallactic_angle) at a latitude.
    ///
    /// Uses ERFA.
    pub fn get_parallactic_angle(self, latitude_rad: f64) -> f64 {
        unsafe { erfa_sys::eraHd2pa(self.ha, self.dec, latitude_rad) }
    }

    /// Get the [parallactic
    /// angle](https://en.wikipedia.org/wiki/Parallactic_angle) at the MWA's
    /// latitude.
    ///
    /// Uses ERFA.
    pub fn get_parallactic_angle_mwa(self) -> f64 {
        unsafe { erfa_sys::eraHd2pa(self.ha, self.dec, MWA_LAT_RAD) }
    }
}

impl std::fmt::Display for HADec {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({}°, {}°)", self.ha.to_degrees(), self.dec.to_degrees())
    }
}

#[cfg(any(test, feature = "approx"))]
impl approx::AbsDiffEq for HADec {
    type Epsilon = f64;

    fn default_epsilon() -> f64 {
        f64::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: f64) -> bool {
        f64::abs_diff_eq(&self.ha, &other.ha, epsilon)
            && f64::abs_diff_eq(&self.dec, &other.dec, epsilon)
    }
}

#[cfg(any(test, feature = "approx"))]
impl approx::RelativeEq for HADec {
    #[inline]
    fn default_max_relative() -> f64 {
        f64::EPSILON
    }

    #[inline]
    fn relative_eq(&self, other: &Self, epsilon: f64, max_relative: f64) -> bool {
        f64::relative_eq(&self.ha, &other.ha, epsilon, max_relative)
            && f64::relative_eq(&self.dec, &other.dec, epsilon, max_relative)
    }

    #[inline]
    fn relative_ne(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        !Self::relative_eq(self, other, epsilon, max_relative)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn to_azel() {
        let hd = HADec::from_degrees(1.0, -35.0);
        let result = hd.to_azel_mwa();
        let expected = AzEl::from_radians(3.240305654530152, 1.425221581624331);
        assert_abs_diff_eq!(result, expected, epsilon = 1e-10);
    }

    #[test]
    fn to_azel2() {
        let hd = HADec::from_degrees(23.0, -35.0);
        let result = hd.to_azel_mwa();
        let expected = AzEl::from_radians(4.215504972991079, 1.1981324538790032);
        assert_abs_diff_eq!(result, expected, epsilon = 1e-10);
    }

    #[test]
    fn separation() {
        let hd1 = HADec::from_degrees(1.0, -35.0);
        let hd2 = HADec::from_degrees(23.0, -35.0);
        let result = hd1.separation(hd2);
        assert_abs_diff_eq!(result, 0.31389018251593337, epsilon = 1e-10);
    }

    #[test]
    fn separation2() {
        let hd1 = HADec::from_degrees(1.0, -35.0);
        let hd2 = HADec::from_degrees(1.1, -35.0);
        let result = hd1.separation(hd2);
        assert_abs_diff_eq!(result, 0.0014296899650293985, epsilon = 1e-10);
    }

    #[test]
    fn separation3() {
        let hd1 = HADec::from_degrees(1.0, -35.0);
        let hd2 = HADec::from_degrees(4.0, 35.0);
        let result = hd1.separation(hd2);
        assert_abs_diff_eq!(result, 1.222708915934097, epsilon = 1e-10);
    }

    #[test]
    fn separation4() {
        let hd1 = HADec::from_degrees(2.0, -35.0);
        let hd2 = HADec::from_degrees(2.0, -35.0);
        let result = hd1.separation(hd2);
        assert_abs_diff_eq!(result, 0.0, epsilon = 1e-10);
    }
}
