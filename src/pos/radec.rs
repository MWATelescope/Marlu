// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Handle (right ascension, declination) coordinates.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use erfa::{
    aliases::eraSeps,
    transform::{cartesian_to_spherical, spherical_to_cartesian},
};

use crate::sexagesimal::{degrees_to_sexagesimal_dms, degrees_to_sexagesimal_hms};

use super::hadec::HADec;
use super::lmn::LMN;

/// A struct containing a Right Ascension and Declination. All units are in
/// radians.
///
/// Note that the serialised units are degrees and are automatically converted
/// when serialising/deserialising.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[allow(clippy::upper_case_acronyms)]
pub struct RADec {
    /// Right ascension \[radians\]
    // TODO: Should RA always be positive?
    #[cfg_attr(feature = "serde", serde(serialize_with = "radians_to_degrees"))]
    #[cfg_attr(feature = "serde", serde(deserialize_with = "degrees_to_radians"))]
    pub ra: f64,

    /// Declination \[radians\]
    #[cfg_attr(feature = "serde", serde(serialize_with = "radians_to_degrees"))]
    #[cfg_attr(feature = "serde", serde(deserialize_with = "degrees_to_radians"))]
    pub dec: f64,
}

#[cfg(feature = "serde")]
fn radians_to_degrees<S: serde::Serializer>(num: &f64, s: S) -> Result<S::Ok, S::Error> {
    s.serialize_f64(num.to_degrees())
}

#[cfg(feature = "serde")]
fn degrees_to_radians<'de, D>(d: D) -> Result<f64, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let num: f64 = serde::Deserialize::deserialize(d)?;
    Ok(num.to_radians())
}

impl RADec {
    /// Make a new [`RADec`] struct from values in radians.
    pub fn from_radians(ra: f64, dec: f64) -> RADec {
        Self { ra, dec }
    }

    /// Make a new [`RADec`] struct from values in degrees.
    pub fn from_degrees(ra: f64, dec: f64) -> RADec {
        Self {
            ra: ra.to_radians(),
            dec: dec.to_radians(),
        }
    }

    /// Make a new [`RADec`] struct from values in radians.
    #[deprecated = "use `RADec::from_radians` instead"]
    pub fn new(ra_rad: f64, dec_rad: f64) -> RADec {
        Self::from_radians(ra_rad, dec_rad)
    }

    /// Make a new [`RADec`] struct from values in degrees.
    #[deprecated = "use `RADec::from_degrees` instead"]
    pub fn new_degrees(ra_deg: f64, dec_deg: f64) -> RADec {
        Self::from_degrees(ra_deg, dec_deg)
    }

    /// Given a local sidereal time, make a new [`HADec`] struct from a [`RADec`].
    pub fn to_hadec(self, lst_rad: f64) -> HADec {
        HADec {
            ha: lst_rad - self.ra,
            dec: self.dec,
        }
    }

    /// Given a local sidereal time, make a new [`RADec`] struct from a [`HADec`].
    pub fn from_hadec(hadec: HADec, lst_rad: f64) -> Self {
        Self {
            ra: lst_rad - hadec.ha,
            dec: hadec.dec,
        }
    }

    /// From a collection of [`RADec`] coordinates and weights, find the average
    /// [`RADec`] position. The lengths of both collection must be the same to
    /// get sensible results. Not providing any [`RADec`] coordinates will make
    /// this function return `None`.
    ///
    /// This function accounts for Right Ascension coordinates that range over
    /// 360 degrees. If the right ascension coordinates span a huge range (e.g.
    /// 0, 90, 180, 270), `None` is returned.
    ///
    /// The average position is found by converting the [`RADec`] coordinates
    /// ("spherical coordinates") to cartesian coordinates (direction cosines),
    /// performing the average, then converting the average back to spherical.
    /// The results have been verified with `scipy.optimize` (a sample program
    /// is provided at the bottom of this file) and online calculators. When
    /// coordinates are widely separated, the average can appear to have strange
    /// values due to them falling on great circles.
    pub fn weighted_average(radecs: &[Self], weights: &[f64]) -> Option<Self> {
        // Accounting for the 360 degree branch cut.
        let mut any_less_than_90 = false;
        let mut any_between_90_270 = false;
        let mut any_greater_than_270 = false;
        for radec in radecs {
            if (0.0..FRAC_PI_2).contains(&radec.ra) {
                any_less_than_90 = true;
            }
            if (FRAC_PI_2..3.0 * FRAC_PI_2).contains(&radec.ra) {
                any_between_90_270 = true;
            }
            if (3.0 * FRAC_PI_2..TAU).contains(&radec.ra) {
                any_greater_than_270 = true;
            }
        }

        let new_cutoff = match (any_less_than_90, any_between_90_270, any_greater_than_270) {
            // All true, this is the Danger Zone; we can't feasibly find the
            // average RADec over a collection of coordinates that span the
            // whole sky.
            //
            // All negative, there are no coordinates.
            (true, true, true) | (false, false, false) => return None,

            // Surrounding 0 or 360.
            (true, false, true) => PI,

            // Easy ones.
            _ => 0.0,
        };

        let mut cart_sum = [0.0; 3];
        let mut weight_sum = 0.0;
        for (&RADec { ra, dec }, w) in radecs.iter().zip(weights.iter()) {
            let ra = if ra > new_cutoff { ra - TAU } else { ra };

            let cart = spherical_to_cartesian(ra, dec);
            cart_sum[0] += cart[0] * w;
            cart_sum[1] += cart[1] * w;
            cart_sum[2] += cart[2] * w;
            weight_sum += w;
        }

        cart_sum[0] /= weight_sum;
        cart_sum[1] /= weight_sum;
        cart_sum[2] /= weight_sum;
        let (ra_average, dec_average) = cartesian_to_spherical(cart_sum);
        let mut average = Self::from_radians(ra_average, dec_average);
        // Keep the RA positive.
        if average.ra < 0.0 {
            average.ra += TAU;
        }

        Some(average)
    }

    /// Get the [LMN] direction cosines from an [`RADec`] and a phase centre.
    ///
    /// Derived using "Coordinate transformations" on page 388 of Synthesis
    /// Imaging in Radio Astronomy II.
    pub fn to_lmn(self, phase_centre: RADec) -> LMN {
        let d_ra = self.ra - phase_centre.ra;
        let (s_d_ra, c_d_ra) = d_ra.sin_cos();
        let (s_dec, c_dec) = self.dec.sin_cos();
        let (pc_s_dec, pc_c_dec) = phase_centre.dec.sin_cos();
        LMN {
            l: c_dec * s_d_ra,
            m: s_dec * pc_c_dec - c_dec * pc_s_dec * c_d_ra,
            n: s_dec * pc_s_dec + c_dec * pc_c_dec * c_d_ra,
        }
    }

    /// Calculate the distance between two sets of coordinates \[radians\].
    pub fn separation(&self, b: Self) -> f64 {
        eraSeps(self.ra, self.dec, b.ra, b.dec)
    }

    /// Given an [`mwalib::MetafitsContext`], make an [`Option<RADec>`] from the
    /// `(ra|dec)_phase_center_degrees` if these are available, otherwise
    /// [`None`].
    #[cfg(feature = "mwalib")]
    pub fn from_mwalib_phase_center(context: &mwalib::MetafitsContext) -> Option<RADec> {
        match (
            context.ra_phase_center_degrees,
            context.dec_phase_center_degrees,
        ) {
            (Some(ra), Some(dec)) => Some(RADec::from_degrees(ra, dec)),
            (..) => None,
        }
    }

    /// Given an [`mwalib::MetafitsContext`], make a [`RADec`] from the
    /// `(ra|dec)_tile_pointing_degrees`.
    #[cfg(feature = "mwalib")]
    pub fn from_mwalib_tile_pointing(context: &mwalib::MetafitsContext) -> RADec {
        RADec::from_degrees(
            context.ra_tile_pointing_degrees,
            context.dec_tile_pointing_degrees,
        )
    }

    /// Given an [`mwalib::MetafitsContext`], make a [`RADec`] from the
    /// `(ra|dec)_phase_center_degrees` if these are available, otherwise use
    /// the `(ra|dec)_tile_pointing_degrees`.
    #[cfg(feature = "mwalib")]
    pub fn from_mwalib_phase_or_pointing(context: &mwalib::MetafitsContext) -> RADec {
        match RADec::from_mwalib_phase_center(context) {
            Some(radec) => radec,
            None => RADec::from_mwalib_tile_pointing(context),
        }
    }
}

impl std::fmt::Display for RADec {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "({:.4}°, {:.4}°) => ({}, {})",
            self.ra.to_degrees(),
            self.dec.to_degrees(),
            degrees_to_sexagesimal_hms(self.ra.to_degrees()),
            degrees_to_sexagesimal_dms(self.dec.to_degrees())
        )
    }
}

#[cfg(any(test, feature = "approx"))]
impl approx::AbsDiffEq for RADec {
    type Epsilon = f64;

    fn default_epsilon() -> f64 {
        f64::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: f64) -> bool {
        f64::abs_diff_eq(&self.ra, &other.ra, epsilon)
            && f64::abs_diff_eq(&self.dec, &other.dec, epsilon)
    }
}

#[cfg(any(test, feature = "approx"))]
impl approx::RelativeEq for RADec {
    #[inline]
    fn default_max_relative() -> f64 {
        f64::EPSILON
    }

    #[inline]
    fn relative_eq(&self, other: &Self, epsilon: f64, max_relative: f64) -> bool {
        f64::relative_eq(&self.ra, &other.ra, epsilon, max_relative)
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
    fn test_to_lmn() {
        let radec = RADec::from_degrees(62.0, -27.5);
        let phase_centre = RADec::from_degrees(60.0, -27.0);
        let lmn = radec.to_lmn(phase_centre);
        let expected = LMN {
            l: 0.03095623164758603,
            m: -0.008971846102111436,
            n: 0.9994804738961642,
        };
        assert_abs_diff_eq!(lmn, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_weighted_pos() {
        // Simple case: both components have a weight of 1.0.
        let radecs = [
            RADec::from_degrees(10.0, 9.0),
            RADec::from_degrees(11.0, 10.0),
        ];
        let w1 = 1.0;
        let w2 = 1.0;
        let result = RADec::weighted_average(&radecs, &[w1, w2]);
        assert!(result.is_some());
        let weighted_pos = result.unwrap();
        assert_abs_diff_eq!(
            weighted_pos,
            RADec {
                ra: 10.49926979_f64.to_radians(),
                dec: 9.500355150_f64.to_radians()
            },
            epsilon = 1e-10
        );

        // Complex case: both components have different weights.
        let w1 = 3.0;
        let result = RADec::weighted_average(&radecs, &[w1, w2]);
        assert!(result.is_some());
        let weighted_pos = result.unwrap();
        assert_abs_diff_eq!(
            weighted_pos,
            RADec {
                ra: 10.24944800_f64.to_radians(),
                dec: 9.250254459_f64.to_radians()
            },
            epsilon = 1e-10
        );
    }

    #[test]
    // This time, make the coordinates go across the 360 degree branch cut.
    fn test_weighted_pos2() {
        let radecs = [
            RADec::from_degrees(60.0, 9.0),
            RADec::from_degrees(359.0, 10.0),
        ];
        let result = RADec::weighted_average(&radecs, &vec![1.0; radecs.len()]);
        assert!(result.is_some());
        let weighted_pos = result.unwrap();
        assert_abs_diff_eq!(
            weighted_pos,
            RADec {
                ra: 29.54928740_f64.to_radians(),
                dec: 10.99095249_f64.to_radians()
            },
            epsilon = 1e-10
        );
    }

    #[test]
    // This time, make the coordinates go from 0 to 260.
    fn test_weighted_pos3() {
        let radecs = [
            RADec::from_degrees(0.0, 10.0),
            RADec::from_degrees(260.0, 10.0),
        ];
        let result = RADec::weighted_average(&radecs, &vec![1.0; radecs.len()]);
        assert!(result.is_some());
        let weighted_pos = result.unwrap();
        assert_abs_diff_eq!(
            weighted_pos,
            RADec {
                ra: 310.0_f64.to_radians(),
                dec: 15.33981450_f64.to_radians()
            },
            epsilon = 1e-10
        );
    }

    #[test]
    // This time, make the coordinates go from 100 to 359.
    fn test_weighted_pos4() {
        let radecs = [
            RADec::from_degrees(100.0, 10.0),
            RADec::from_degrees(359.0, 10.0),
        ];
        let result = RADec::weighted_average(&radecs, &vec![1.0; radecs.len()]);
        assert!(result.is_some());
        let weighted_pos = result.unwrap();
        assert_abs_diff_eq!(
            weighted_pos,
            RADec {
                ra: 49.5_f64.to_radians(),
                dec: 15.49388365_f64.to_radians()
            },
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_weighted_pos5() {
        let radecs = [
            RADec::from_degrees(0.0, 0.0),
            RADec::from_degrees(90.0, 0.0),
            RADec::from_degrees(180.0, 0.0),
            RADec::from_degrees(270.0, 0.0),
        ];
        let result = RADec::weighted_average(&radecs, &vec![1.0; radecs.len()]);
        assert!(result.is_none());
    }

    #[test]
    fn test_weighted_pos_single() {
        let c = RADec::from_radians(0.5, 0.75);
        let w = 1.0;
        let result = RADec::weighted_average(&[c], &[w]);
        assert!(result.is_some());
        let weighted_pos = result.unwrap();
        assert_abs_diff_eq!(weighted_pos, RADec { ra: 0.5, dec: 0.75 }, epsilon = 1e-10);
    }

    #[test]
    fn test_weighted_pos_empty() {
        let arr: Vec<RADec> = vec![];
        assert!(RADec::weighted_average(&arr, &[1.0]).is_none());
    }

    #[test]
    fn test_weighted_pos_all_sky() {
        let radecs = [
            RADec::from_degrees(10.0, 0.0),
            RADec::from_degrees(80.0, 0.0),
            RADec::from_degrees(150.0, 0.0),
            RADec::from_degrees(280.0, 0.0),
        ];
        assert!(RADec::weighted_average(&radecs, &vec![1.0; radecs.len()]).is_none());
    }

    #[test]
    fn test_display_radec() {
        let radec = RADec { ra: 0.0, dec: 0.0 };
        let result = format!("{radec}");
        assert!(!result.is_empty());
    }

    #[test]
    #[cfg(feature = "serde")]
    fn test_serde() {
        let radec = RADec::from_degrees(60.0, -30.0);
        let result = serde_json::to_string(&radec);
        assert!(result.is_ok(), "{:?}", result.err());
        let json = result.unwrap();

        let result = serde_json::from_str(&json);
        assert!(result.is_ok(), "{:?}", result.err());
        let radec2 = result.unwrap();

        assert_abs_diff_eq!(radec, radec2);

        // No issue with the pretty renderer.
        let result = serde_json::to_string_pretty(&radec);
        assert!(result.is_ok(), "{:?}", result.err());
        let json = result.unwrap();

        let result = serde_json::from_str(&json);
        assert!(result.is_ok(), "{:?}", result.err());
        let radec2 = result.unwrap();

        assert_abs_diff_eq!(radec, radec2);
    }

    #[test]
    #[cfg(feature = "serde")]
    fn test_deserialise_json() {
        let json = "{\"ra\": 1.23, \"dec\": -0.57}";

        let result = serde_json::from_str(&json);
        assert!(result.is_ok(), "{:?}", result.err());
        let radec: RADec = result.unwrap();

        assert_abs_diff_eq!(
            radec,
            RADec {
                ra: 0.021467549799530253,
                dec: -0.009948376736367677
            }
        );
    }
}

/* Sample Python program to find the average RADec from a collection of RADecs.

#!/usr/bin/env python

import numpy as np
from scipy.optimize import minimize
from astropy.coordinates import SkyCoord
from astropy import units as u

# radecs = [(0.0, 0.0), (1.0, 0.0), (0.5, 0.5)]
# guess = (0.5, 0.25)

# radecs = [(0.0, 45.0), (10.0, 45.0)]
# guess = (0.0, 75.0)

# radecs = [(6.454873596315, -26.036398659571013),
#           (6.45640148205, -26.043503854016866),
#           (6.4706449552650005, -26.020419803182456),
#           (6.4869815745, -26.04649083746405)
#           ]
# guess = (6.45, -26.0)

radecs = [(0, 10), (10, 10)]
guess = (np.mean([radec[0] for radec in radecs]), np.mean([radec[1] for radec in radecs]))

print("radecs:", radecs)
print("guess:", guess)

def objective(x0):
    cost = 0.0
    average = SkyCoord(ra=x0[0]*u.deg, dec=x0[1]*u.deg)
    for radec in radecs:
        radec = SkyCoord(ra=radec[0]*u.deg, dec=radec[1]*u.deg)
        sep = radec.separation(average).value**2
        cost += sep
    return cost

result = minimize(objective, guess, tol=1e-10)
# print(result)
print(result.x)

 */
