// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Handle East, North and Height coordinates (typically associated with MWA
//! tiles).

use crate::{constants::MWA_LAT_RAD, XyzGeodetic};

/// East, North and Height coordinates.
#[derive(Clone, Copy, Debug, Default)]
#[allow(clippy::upper_case_acronyms)]
pub struct ENH {
    /// East \[metres\]
    pub e: f64,
    /// North \[metres\]
    pub n: f64,
    /// Height \[metres\]
    pub h: f64,
}

impl ENH {
    /// Convert coords in local topocentric East, North, Height units to 'local'
    /// [`XyzGeodetic`] units. Local means Z points north, X points through the
    /// equator from the geocenter along the local meridian and Y is East. This
    /// is like the absolute system except that zero longitude is now the local
    /// meridian rather than prime meridian. Latitude is geodetic, in radians.
    /// This is what you want for constructing the local antenna positions in a
    /// UVFITS antenna table.
    ///
    /// Taken from the third edition of Interferometry and Synthesis in Radio
    /// Astronomy, chapter 4: Geometrical Relationships, Polarimetry, and the
    /// Measurement Equation.
    pub fn to_xyz(self, latitude_rad: f64) -> XyzGeodetic {
        let (s_lat, c_lat) = latitude_rad.sin_cos();
        Self::to_xyz_inner(self, s_lat, c_lat)
    }

    /// Convert coords in local topocentric East, North, Height units to 'local'
    /// [`XyzGeodetic`] units. See [`ENH::to_xyz`] for more information. This
    /// function is less convenient than [`ENH::to_xyz`], but is slightly more
    /// efficient because the caller can prevent needless `sin` and `cos`
    /// calculations.
    pub fn to_xyz_inner(self, sin_latitude: f64, cos_latitude: f64) -> XyzGeodetic {
        XyzGeodetic {
            x: -self.n * sin_latitude + self.h * cos_latitude,
            y: self.e,
            z: self.n * cos_latitude + self.h * sin_latitude,
        }
    }

    /// Convert [`ENH`] coordinates to [`XyzGeodetic`] for the MWA's latitude.
    pub fn to_xyz_mwa(self) -> XyzGeodetic {
        self.to_xyz(MWA_LAT_RAD)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn convert_enh_to_xyz_test() {
        let enh = ENH {
            n: -101.530,
            e: -585.675,
            h: 375.212,
        };
        let xyz = ENH::to_xyz_mwa(enh);
        assert_abs_diff_eq!(xyz.x, 289.56928486613185, epsilon = 1e-10);
        assert_abs_diff_eq!(xyz.y, -585.675, epsilon = 1e-10);
        assert_abs_diff_eq!(xyz.z, -259.3106536687549, epsilon = 1e-10);
    }
}
