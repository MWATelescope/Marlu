// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Handle (azimuth, elevation) coordinates (also known as horizontal
//! coordinates).

use std::f64::consts::FRAC_PI_2;
use super::hadec::HADec;

/// A struct containing an Azimuth and Elevation. All units are in radians.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AzEl {
    /// Hour angle \[radians\]
    pub az: f64,
    /// Declination \[radians\]
    pub el: f64,
}

impl AzEl {
    /// Make a new [AzEl] struct from values in radians.
    pub fn new(az_rad: f64, el_rad: f64) -> AzEl {
        Self {
            az: az_rad,
            el: el_rad,
        }
    }

    /// Make a new [AzEl] struct from values in degrees.
    pub fn new_degrees(az_deg: f64, el_deg: f64) -> AzEl {
        Self::new(az_deg.to_radians(), el_deg.to_radians())
    }

    /// Get the zenith angle in radians.
    pub fn za(self) -> f64 {
        FRAC_PI_2 - self.el
    }

    /// Convert the horizon coordinates to equatorial coordinates (Hour Angle
    /// and Declination), given the local latitude on Earth.
    ///
    /// Uses ERFA.
    pub fn to_hadec(self, latitude_rad: f64) -> HADec {
        let mut ha = 0.0;
        let mut dec = 0.0;
        unsafe { erfa_sys::eraAe2hd(self.az, self.el, latitude_rad, &mut ha, &mut dec) }
        HADec::new(ha, dec)
    }

    /// Convert the horizon coordinates to equatorial coordinates (Hour Angle
    /// and Declination) for the MWA's location.
    ///
    /// Uses ERFA.
    pub fn to_hadec_mwa(self) -> HADec {
        self.to_hadec(crate::constants::MWA_LAT_RAD)
    }
}

impl std::fmt::Display for AzEl {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "({:.4}°, {:.4}°)",
            self.az.to_degrees(),
            self.el.to_degrees()
        )
    }
}

#[cfg(test)]
use crate::approx::AbsDiffEq;

#[cfg(test)]
impl AbsDiffEq for AzEl {
    type Epsilon = f64;

    fn default_epsilon() -> f64 {
        f64::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: f64) -> bool {
        f64::abs_diff_eq(&self.az, &other.az, epsilon)
            && f64::abs_diff_eq(&self.el, &other.el, epsilon)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::approx::assert_abs_diff_eq;

    #[test]
    fn to_hadec() {
        let ae = AzEl::new_degrees(45.0, 30.0);
        let result = ae.to_hadec(-0.497600);
        let expected = HADec::new(-0.6968754873551053, 0.3041176697804004);
        assert_abs_diff_eq!(result, expected, epsilon = 1e-10);
    }

    #[test]
    fn to_hadec2() {
        let ae = AzEl::new(0.261700, 0.785400);
        let result = ae.to_hadec(-0.897600);
        let expected = HADec::new(-0.185499449332533, -0.12732312479328656);
        assert_abs_diff_eq!(result, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_za() {
        let ae = AzEl::new(0.261700, 0.785400);
        let za = ae.za();
        assert_abs_diff_eq!(za, 0.7853963268, epsilon = 1e-10);
    }
}
