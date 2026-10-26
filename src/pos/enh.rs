// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Handle East, North and Height coordinates (typically associated with MWA
//! tiles).

#[cfg(feature = "mwalib")]
use crate::constants::MWA_HEIGHT_M;
use crate::{constants::MWA_LAT_RAD, LatLngHeight, XyzGeodetic};

/// East, North and Height coordinates.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
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
    /// The position of a metafits antenna relative to the MWA's array centre.
    /// Metafits heights are altitudes above sea level, so [`MWA_HEIGHT_M`] is
    /// taken off them; using them as they are would count the array's height
    /// twice when converting to geocentric coordinates (as `cotter` did,
    /// putting every tile ~378 m too high).
    ///
    /// The MWA database (and so the metafits) stores the 128T and hex tiles'
    /// heights as their local "up" plus the array's altitude, but the
    /// long-baseline tiles (`LB*`) as their raw surveyed heights above the
    /// ellipsoid. The latter differ from "up" by the Earth's curvature (~1 m
    /// at 3.5 km), so they are converted with
    /// [`ENH::from_ellipsoidal_height`]. See the MWA wiki page "MWA Tile
    /// coordinate errors".
    #[cfg(feature = "mwalib")]
    pub fn from_mwalib_antenna(antenna: &mwalib::Antenna) -> ENH {
        if antenna.tile_name.starts_with("LB") {
            ENH::from_ellipsoidal_height(
                antenna.east_m,
                antenna.north_m,
                antenna.height_m,
                LatLngHeight::mwa(),
            )
        } else {
            ENH {
                e: antenna.east_m,
                n: antenna.north_m,
                h: antenna.height_m - MWA_HEIGHT_M,
            }
        }
    }

    /// The position of a point with local East and North coordinates relative
    /// to `array_pos`, and a height above the WGS84 ellipsoid. The returned
    /// height is along `array_pos`'s vertical, i.e. it includes the drop of
    /// the Earth's curvature (approximately (e² + n²) / 2R) below the array's
    /// tangent plane.
    pub fn from_ellipsoidal_height(
        e: f64,
        n: f64,
        ellipsoidal_height_m: f64,
        array_pos: LatLngHeight,
    ) -> ENH {
        let (s_lat, c_lat) = array_pos.latitude_rad.sin_cos();
        let mut enh = ENH {
            e,
            n,
            h: ellipsoidal_height_m - array_pos.height_metres,
        };
        // d(ellipsoidal height)/dh is ~1, so this converges to well under a
        // micrometre within a few iterations.
        for _ in 0..4 {
            let height = enh
                .to_xyz_inner(s_lat, c_lat)
                .to_geocentric(array_pos)
                .to_earth_wgs84()
                .height_metres;
            enh.h += ellipsoidal_height_m - height;
        }
        enh
    }

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

#[cfg(any(test, feature = "approx"))]
impl approx::AbsDiffEq for ENH {
    type Epsilon = f64;

    fn default_epsilon() -> f64 {
        f64::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: f64) -> bool {
        f64::abs_diff_eq(&self.e, &other.e, epsilon)
            && f64::abs_diff_eq(&self.n, &other.n, epsilon)
            && f64::abs_diff_eq(&self.h, &other.h, epsilon)
    }
}

#[cfg(any(test, feature = "approx"))]
impl approx::RelativeEq for ENH {
    #[inline]
    fn default_max_relative() -> f64 {
        f64::EPSILON
    }

    #[inline]
    fn relative_eq(&self, other: &Self, epsilon: f64, max_relative: f64) -> bool {
        f64::relative_eq(&self.e, &other.e, epsilon, max_relative)
            && f64::relative_eq(&self.n, &other.n, epsilon, max_relative)
            && f64::relative_eq(&self.h, &other.h, epsilon, max_relative)
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
    fn convert_enh_to_xyz_test() {
        let enh = ENH {
            n: -101.530,
            e: -585.675,
            h: 375.212,
        };
        let xyz = ENH::to_xyz_mwa(enh);
        assert_abs_diff_eq!(
            xyz,
            XyzGeodetic {
                x: 289.56928486613185,
                y: -585.675,
                z: -259.3106536687549
            },
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_from_ellipsoidal_height() {
        // At the array centre, nothing changes.
        let array = LatLngHeight::mwa();
        let enh = ENH::from_ellipsoidal_height(0.0, 0.0, array.height_metres, array);
        assert_abs_diff_eq!(enh.h, 0.0, epsilon = 1e-6);

        // LBD7, LBA1 and LBB2, with their raw survey heights; the expected "up"
        // values (survey heights converted to local cartesian coordinates)
        // are from the MWA wiki page "MWA Tile coordinate errors", rounded
        // to the nearest millimetre.
        for (e, n, height, up) in [
            (3105.64, 1640.33, 376.795, 375.828),
            (-1999.81, 206.85, 372.921, 372.604),
            (-330.61, 2829.62, 373.367, 372.728),
        ] {
            let enh = ENH::from_ellipsoidal_height(e, n, height, LatLngHeight::mwa());
            assert_abs_diff_eq!(enh.e, e);
            assert_abs_diff_eq!(enh.n, n);
            assert_abs_diff_eq!(enh.h, up - 377.8269, epsilon = 2e-3);
        }
    }

    #[test]
    #[cfg(feature = "mwalib")]
    fn test_from_mwalib_antenna_long_baseline() {
        let context =
            mwalib::MetafitsContext::new("tests/data/1196175296_mwa_ord/1196175296.metafits", None)
                .unwrap();
        let lbd7 = context
            .antennas
            .iter()
            .find(|a| a.tile_name == "LBD7")
            .unwrap();
        // The metafits has the raw survey height (376.80 m) for this tile; its
        // local "up" is ~0.97 m lower than that minus the array's height.
        let enh = ENH::from_mwalib_antenna(lbd7);
        assert_abs_diff_eq!(enh.h, 375.828 - 377.8269, epsilon = 1e-2);
        assert!(lbd7.height_m - MWA_HEIGHT_M - enh.h > 0.9);
    }
}
