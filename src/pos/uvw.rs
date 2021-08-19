// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Handle UVW coordinates.

use super::hadec::HADec;
use super::xyz::XyzGeodetic;

/// The (u,v,w) coordinates of a baseline. All units are in terms of wavelength,
/// with units of metres.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[allow(clippy::upper_case_acronyms)]
pub struct UVW {
    /// u coordinate \[meters\]
    pub u: f64,
    /// v coordinate \[meters\]
    pub v: f64,
    /// w coordinate \[meters\]
    pub w: f64,
}

impl UVW {
    /// Convert an [XyzGeodetic] to [UVW], given the phase centre.
    ///
    /// This is Equation 4.1 of: Interferometry and Synthesis in Radio
    /// Astronomy, Third Edition, Section 4: Geometrical Relationships,
    /// Polarimetry, and the Measurement Equation.
    pub fn from_xyz(xyz: XyzGeodetic, phase_centre: HADec) -> Self {
        let (s_ha, c_ha) = phase_centre.ha.sin_cos();
        let (s_dec, c_dec) = phase_centre.dec.sin_cos();
        Self::from_xyz_inner(xyz, s_ha, c_ha, s_dec, c_dec)
    }

    /// Convert an [XyzGeodetic] to [UVW], given the phase centre. This function
    /// is less convenient than [UVW::from_xyz()], but may be better in tight
    /// loops as the `sin` and `cos` of the phase centre don't need to be
    /// uselessly re-calculated.
    ///
    /// This is Equation 4.1 of: Interferometry and Synthesis in Radio
    /// Astronomy, Third Edition, Section 4: Geometrical Relationships,
    /// Polarimetry, and the Measurement Equation.
    pub fn from_xyz_inner(xyz: XyzGeodetic, s_ha: f64, c_ha: f64, s_dec: f64, c_dec: f64) -> Self {
        Self {
            u: s_ha * xyz.x + c_ha * xyz.y,
            v: -s_dec * c_ha * xyz.x + s_dec * s_ha * xyz.y + c_dec * xyz.z,
            w: c_dec * c_ha * xyz.x - c_dec * s_ha * xyz.y + s_dec * xyz.z,
        }
    }
}

impl std::ops::Sub<UVW> for UVW {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        UVW {
            u: self.u - rhs.u,
            v: self.v - rhs.v,
            w: self.w - rhs.w,
        }
    }
}

impl std::ops::Mul<f64> for UVW {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self {
        UVW {
            u: self.u * rhs,
            v: self.v * rhs,
            w: self.w * rhs,
        }
    }
}

impl std::ops::Div<f64> for UVW {
    type Output = Self;

    fn div(self, rhs: f64) -> Self {
        UVW {
            u: self.u / rhs,
            v: self.v / rhs,
            w: self.w / rhs,
        }
    }
}

#[cfg(test)]
impl approx::AbsDiffEq for UVW {
    type Epsilon = f64;

    fn default_epsilon() -> f64 {
        f64::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: f64) -> bool {
        f64::abs_diff_eq(&self.u, &other.u, epsilon)
            && f64::abs_diff_eq(&self.v, &other.v, epsilon)
            && f64::abs_diff_eq(&self.w, &other.w, epsilon)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::*;

    #[test]
    fn test_uvw_mul() {
        let uvw = UVW {
            u: 1.0,
            v: 2.0,
            w: 3.0,
        } * 3.0;
        assert_abs_diff_eq!(uvw.u, 3.0);
        assert_abs_diff_eq!(uvw.v, 6.0);
        assert_abs_diff_eq!(uvw.w, 9.0);
    }

    #[test]
    fn test_uvw_div() {
        let uvw = UVW {
            u: 3.0,
            v: 6.0,
            w: 9.0,
        } / 3.0;
        assert_abs_diff_eq!(uvw.u, 1.0);
        assert_abs_diff_eq!(uvw.v, 2.0);
        assert_abs_diff_eq!(uvw.w, 3.0);
    }
}
