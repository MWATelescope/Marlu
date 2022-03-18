// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Handle (x,y,z) coordinates of an antenna (a.k.a. tile or station), geodetic
//! or geocentric.
//!
//! hyperdrive prefers to keep track of [XyzGeodetic] coordinates, as these are
//! what are needed to calculate [UVW]s.
//!
//! This coordinate system is discussed at length in Interferometry and
//! Synthesis in Radio Astronomy, Third Edition, Section 4: Geometrical
//! Relationships, Polarimetry, and the Measurement Equation.

// TODO: Account for northing and eastings. Australia drifts by ~7cm/year, and
// the ellipsoid model probably need to be changed too!

use rayon::prelude::*;

use super::ErfaError;
use crate::{
    constants::MWA_LAT_RAD,
    math::{baseline_to_tiles, cross_correlation_baseline_to_tiles},
    HADec, LatLngHeight, ENH, UVW,
};

/// The geodetic (x,y,z) coordinates of an antenna (a.k.a. tile or station). All
/// units are in metres.
///
/// This coordinate system is discussed at length in Interferometry and
/// Synthesis in Radio Astronomy, Third Edition, Section 4: Geometrical
/// Relationships, Polarimetry, and the Measurement Equation.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct XyzGeodetic {
    /// x-coordinate \[meters\]
    pub x: f64,
    /// y-coordinate \[meters\]
    pub y: f64,
    /// z-coordinate \[meters\]
    pub z: f64,
}

impl XyzGeodetic {
    /// Convert [XyzGeodetic] coordinates at a latitude to [ENH] coordinates.
    pub fn to_enh(self, latitude: f64) -> ENH {
        let (s_lat, c_lat) = latitude.sin_cos();
        Self::to_enh_inner(self, s_lat, c_lat)
    }

    /// Convert [XyzGeodetic] coordinates at a latitude to [ENH] coordinates.
    /// This function is less convenient than [XyzGeodetic::to_enh()], but is
    /// slightly more efficient because the caller can prevent needless `sin`
    /// and `cos` calculations.
    pub fn to_enh_inner(self, sin_latitude: f64, cos_latitude: f64) -> ENH {
        ENH {
            e: self.y,
            n: -self.x * sin_latitude + self.z * cos_latitude,
            h: self.x * cos_latitude + self.z * sin_latitude,
        }
    }

    /// Convert [XyzGeodetic] coordinates at the MWA's latitude to [ENH]
    /// coordinates.
    pub fn to_enh_mwa(self) -> ENH {
        self.to_enh(MWA_LAT_RAD)
    }

    /// Convert a [XyzGeodetic] coordinate to [XyzGeocentric].
    pub fn to_geocentric(self, earth_pos: LatLngHeight) -> Result<XyzGeocentric, ErfaError> {
        let (sin_longitude, cos_longitude) = earth_pos.longitude_rad.sin_cos();
        let geocentric_vector = XyzGeocentric::get_geocentric_vector(earth_pos)?;
        Ok(XyzGeodetic::to_geocentric_inner(
            self,
            geocentric_vector,
            sin_longitude,
            cos_longitude,
        ))
    }

    /// Convert a [XyzGeodetic] coordinate to [XyzGeocentric]. This function is
    /// less convenient than [XyzGeodetic::to_geocentric], but may be better in
    /// tight loops as the arguments to this function don't need to be uselessly
    /// re-calculated.
    pub fn to_geocentric_inner(
        self,
        geocentric_vector: XyzGeocentric,
        sin_longitude: f64,
        cos_longitude: f64,
    ) -> XyzGeocentric {
        let xtemp = self.x * cos_longitude - self.y * sin_longitude;
        let y = self.x * sin_longitude + self.y * cos_longitude;
        let x = xtemp;

        XyzGeocentric {
            x: x + geocentric_vector.x,
            y: y + geocentric_vector.y,
            z: self.z + geocentric_vector.z,
        }
    }

    /// Convert a [XyzGeodetic] coordinate to [XyzGeocentric], using the MWA's
    /// location.
    pub fn to_geocentric_mwa(self) -> Result<XyzGeocentric, ErfaError> {
        self.to_geocentric(LatLngHeight::new_mwa())
    }

    /// For each tile listed in an [`mwalib::MetafitsContext`], calculate a
    /// [XyzGeodetic] coordinate.
    ///
    /// Note that the RF inputs are ordered by antenna number, **not** the
    /// "input"; e.g. in the metafits file, Tile104 is often the first tile
    /// listed ("input" 0), Tile103 second ("input" 2), so the first baseline
    /// would naively be between Tile104 and Tile103.
    #[cfg(feature = "mwalib")]
    pub fn get_tiles(context: &mwalib::MetafitsContext, latitude_rad: f64) -> Vec<XyzGeodetic> {
        let (sin_lat, cos_lat) = latitude_rad.sin_cos();
        context
            .rf_inputs
            .iter()
            // There is an RF input for both tile polarisations. The ENH
            // coordinates are the same for both polarisations of a tile; ignore
            // the RF input if it's associated with Y.
            .filter(|rf| matches!(rf.pol, mwalib::Pol::Y))
            .map(|rf| {
                ENH {
                    e: rf.east_m,
                    n: rf.north_m,
                    h: rf.height_m,
                }
                .to_xyz_inner(sin_lat, cos_lat)
            })
            .collect()
    }

    /// For each tile listed in an [`mwalib::MetafitsContext`], calculate a
    /// [XyzGeodetic] coordinate assuming the MWA's latitude.
    ///
    /// Note that the RF inputs are ordered by antenna number, **not** the
    /// "input"; e.g. in the metafits file, Tile104 is often the first tile
    /// listed ("input" 0), Tile103 second ("input" 2), so the first baseline
    /// would naively be between Tile104 and Tile103.
    #[cfg(feature = "mwalib")]
    pub fn get_tiles_mwa(context: &mwalib::MetafitsContext) -> Vec<XyzGeodetic> {
        Self::get_tiles(context, MWA_LAT_RAD)
    }
}

/// Convert [XyzGeodetic] tile coordinates to [UVW] baseline coordinates without
/// having to form [XyzGeodetic] baselines first.
pub fn xyzs_to_uvws(xyzs: &[XyzGeodetic], phase_centre: HADec) -> Vec<UVW> {
    let (s_ha, c_ha) = phase_centre.ha.sin_cos();
    let (s_dec, c_dec) = phase_centre.dec.sin_cos();
    // Get a UVW for each tile.
    let tile_uvws: Vec<UVW> = xyzs
        .iter()
        .map(|&xyz| UVW::from_xyz_inner(xyz, s_ha, c_ha, s_dec, c_dec))
        .collect();
    // Take the difference of every pair of UVWs.
    let num_tiles = xyzs.len();
    let num_baselines = (num_tiles * (num_tiles + 1)) / 2;
    let mut bl_uvws = Vec::with_capacity(num_baselines);
    for i in 0..num_tiles {
        for j in i..num_tiles {
            bl_uvws.push(tile_uvws[i] - tile_uvws[j]);
        }
    }
    bl_uvws
}

/// Convert [XyzGeodetic] tile coordinates to [UVW] baseline coordinates without
/// having to form [XyzGeodetic] baselines first. This function performs
/// calculations in parallel.
pub fn xyzs_to_uvws_parallel(xyzs: &[XyzGeodetic], phase_centre: HADec) -> Vec<UVW> {
    let (s_ha, c_ha) = phase_centre.ha.sin_cos();
    let (s_dec, c_dec) = phase_centre.dec.sin_cos();
    // Get a UVW for each tile.
    let tile_uvws: Vec<UVW> = xyzs
        .par_iter()
        .map(|&xyz| UVW::from_xyz_inner(xyz, s_ha, c_ha, s_dec, c_dec))
        .collect();
    // Take the difference of every pair of UVWs.
    let num_tiles = xyzs.len();
    let num_baselines = (num_tiles * (num_tiles + 1)) / 2;
    (0..num_baselines)
        .into_par_iter()
        .map(|i_bl| {
            let (i, j) = baseline_to_tiles(num_tiles, i_bl);
            tile_uvws[i] - tile_uvws[j]
        })
        .collect()
}

/// Convert [XyzGeodetic] tile coordinates to [UVW] baseline coordinates without
/// having to form [XyzGeodetic] baselines first. Cross-correlation baselines
/// only.
pub fn xyzs_to_cross_uvws(xyzs: &[XyzGeodetic], phase_centre: HADec) -> Vec<UVW> {
    let (s_ha, c_ha) = phase_centre.ha.sin_cos();
    let (s_dec, c_dec) = phase_centre.dec.sin_cos();
    // Get a UVW for each tile.
    let tile_uvws: Vec<UVW> = xyzs
        .iter()
        .map(|&xyz| UVW::from_xyz_inner(xyz, s_ha, c_ha, s_dec, c_dec))
        .collect();
    // Take the difference of every pair of UVWs.
    let num_tiles = xyzs.len();
    let num_baselines = (num_tiles * (num_tiles - 1)) / 2;
    let mut bl_uvws = Vec::with_capacity(num_baselines);
    for i in 0..num_tiles {
        for j in i + 1..num_tiles {
            bl_uvws.push(tile_uvws[i] - tile_uvws[j]);
        }
    }
    bl_uvws
}

/// Convert [XyzGeodetic] tile coordinates to [UVW] baseline coordinates without
/// having to form [XyzGeodetic] baselines first. This function performs
/// calculations in parallel. Cross-correlation baselines only.
pub fn xyzs_to_cross_uvws_parallel(xyzs: &[XyzGeodetic], phase_centre: HADec) -> Vec<UVW> {
    let (s_ha, c_ha) = phase_centre.ha.sin_cos();
    let (s_dec, c_dec) = phase_centre.dec.sin_cos();
    // Get a UVW for each tile.
    let tile_uvws: Vec<UVW> = xyzs
        .par_iter()
        .map(|&xyz| UVW::from_xyz_inner(xyz, s_ha, c_ha, s_dec, c_dec))
        .collect();
    // Take the difference of every pair of UVWs.
    let num_tiles = xyzs.len();
    let num_baselines = (num_tiles * (num_tiles - 1)) / 2;
    (0..num_baselines)
        .into_par_iter()
        .map(|i_bl| {
            let (i, j) = cross_correlation_baseline_to_tiles(num_tiles, i_bl);
            tile_uvws[i] - tile_uvws[j]
        })
        .collect()
}

impl std::ops::Sub<XyzGeodetic> for XyzGeodetic {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        XyzGeodetic {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

/// The geocentric (x,y,z) coordinates of an antenna (a.k.a. tile or station).
/// All units are in metres.
///
/// This coordinate system is discussed at length in Interferometry and
/// Synthesis in Radio Astronomy, Third Edition, Section 4: Geometrical
/// Relationships, Polarimetry, and the Measurement Equation.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct XyzGeocentric {
    /// x-coordinate \[meters\]
    pub x: f64,
    /// y-coordinate \[meters\]
    pub y: f64,
    /// z-coordinate \[meters\]
    pub z: f64,
}

impl XyzGeocentric {
    /// Get a geocentric coordinate vector with the given geodetic coordinates
    /// (longitude, latitude and height). The ellipsoid model is WGS84.
    pub fn get_geocentric_vector(earth_pos: LatLngHeight) -> Result<XyzGeocentric, ErfaError> {
        let mut geocentric_vector: [f64; 3] = [0.0; 3];
        let status = unsafe {
            erfa_sys::eraGd2gc(
                erfa_sys::ERFA_WGS84 as i32,    // ellipsoid identifier (Note 1)
                earth_pos.longitude_rad,        // longitude (radians, east +ve)
                earth_pos.latitude_rad,         // latitude (geodetic, radians, Note 3)
                earth_pos.height_metres,        // height above ellipsoid (geodetic, Notes 2,3)
                geocentric_vector.as_mut_ptr(), // geocentric vector (Note 2)
            )
        };
        if status != 0 {
            return Err(ErfaError {
                source_file: file!(),
                source_line: line!(),
                status,
                function: "eraGd2gc",
            });
        }
        Ok(XyzGeocentric {
            x: geocentric_vector[0],
            y: geocentric_vector[1],
            z: geocentric_vector[2],
        })
    }

    /// Get a geocentric coordinate vector with the MWA's location. This
    /// function just calls [XyzGeocentric::get_geocentric_vector] with
    /// [MWA_LONG_RAD], [MWA_LAT_RAD] and [MWA_HEIGHT_M].
    pub fn get_geocentric_vector_mwa() -> Result<XyzGeocentric, ErfaError> {
        Self::get_geocentric_vector(LatLngHeight::new_mwa())
    }

    /// Convert a [XyzGeocentric] coordinate to [XyzGeodetic].
    pub fn to_geodetic(self, earth_pos: LatLngHeight) -> Result<XyzGeodetic, ErfaError> {
        let geocentric_vector = XyzGeocentric::get_geocentric_vector(earth_pos)?;
        let (sin_longitude, cos_longitude) = earth_pos.longitude_rad.sin_cos();
        let geodetic =
            XyzGeocentric::to_geodetic_inner(self, geocentric_vector, sin_longitude, cos_longitude);
        Ok(geodetic)
    }

    /// Convert a [XyzGeocentric] coordinate to [XyzGeodetic]. This function is
    /// less convenient than [XyzGeocentric::to_geodetic()], but may be better
    /// in tight loops as the arguments to this function don't need to be
    /// uselessly re-calculated.
    pub fn to_geodetic_inner(
        self,
        geocentric_vector: XyzGeocentric,
        sin_longitude: f64,
        cos_longitude: f64,
    ) -> XyzGeodetic {
        let geodetic = XyzGeodetic {
            x: self.x - geocentric_vector.x,
            y: self.y - geocentric_vector.y,
            z: self.z - geocentric_vector.z,
        };

        let xtemp = geodetic.x * cos_longitude - geodetic.y * -sin_longitude;
        let y = geodetic.x * -sin_longitude + geodetic.y * cos_longitude;
        let x = xtemp;
        XyzGeodetic {
            x,
            y,
            z: geodetic.z,
        }
    }

    /// Convert a [XyzGeocentric] coordinate to [XyzGeodetic], using the MWA's
    /// location.
    pub fn to_geodetic_mwa(self) -> Result<XyzGeodetic, ErfaError> {
        self.to_geodetic(LatLngHeight::new_mwa())
    }
}

#[cfg(test)]
use approx::AbsDiffEq;

#[cfg(test)]
impl AbsDiffEq for XyzGeodetic {
    type Epsilon = f64;

    fn default_epsilon() -> f64 {
        f64::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: f64) -> bool {
        f64::abs_diff_eq(&self.x, &other.x, epsilon)
            && f64::abs_diff_eq(&self.y, &other.y, epsilon)
            && f64::abs_diff_eq(&self.z, &other.z, epsilon)
    }
}

#[cfg(test)]
impl AbsDiffEq for XyzGeocentric {
    type Epsilon = f64;

    fn default_epsilon() -> f64 {
        f64::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: f64) -> bool {
        f64::abs_diff_eq(&self.x, &other.x, epsilon)
            && f64::abs_diff_eq(&self.y, &other.y, epsilon)
            && f64::abs_diff_eq(&self.z, &other.z, epsilon)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ndarray::Array1;
    use approx::assert_abs_diff_eq;

    use crate::constants::{
        COTTER_MWA_HEIGHT_METRES, COTTER_MWA_LATITUDE_RADIANS, COTTER_MWA_LONGITUDE_RADIANS,
    };

    #[test]
    fn test_geocentric_to_geodetic() {
        // Do everything manually.
        let geocentric_vector = XyzGeocentric {
            x: -2559453.2905955315,
            y: 5095371.7354411585,
            z: -2849056.7735717744,
        };
        let sin_longitude = 0.8936001831599957;
        let cos_longitude = -0.44886380190033387;
        let geocentric = XyzGeocentric {
            x: -2559043.7415729975,
            y: 5095823.023550426,
            z: -2849455.5775171486,
        };
        let result = geocentric.to_geodetic_inner(geocentric_vector, sin_longitude, cos_longitude);
        let expected = XyzGeodetic {
            x: 219.43940577989025,
            y: -568.5399780273752,
            z: -398.80394537420943,
        };
        assert_abs_diff_eq!(result, expected);

        // Do everything automatically.
        let result = geocentric
            .to_geodetic(LatLngHeight {
                longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
                latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
                height_metres: COTTER_MWA_HEIGHT_METRES,
            })
            .unwrap();
        assert_abs_diff_eq!(result, expected);
    }

    #[test]
    fn test_geocentric_to_geodetic_and_back() {
        // These geodetic XYZ positions are taken from a uvfits made from cotter
        // for Tile011.
        let uvfits_xyz = XyzGeodetic {
            x: 4.56250049e+02,
            y: -1.49785004e+02,
            z: 6.80459899e+01,
        };
        // These geocentric XYZ positions are taken from a MS made from cotter
        // for Tile011.
        let ms_xyz = XyzGeocentric {
            x: -2559524.23682043,
            y: 5095846.67363471,
            z: -2848988.72758185,
        };

        // Check the conversion of geocentric to geodetic.
        let result = ms_xyz.to_geodetic_mwa();
        assert!(result.is_ok());
        let local_xyz = result.unwrap();

        // cotter's MWA coordinates are a little off of what is in mwalib.
        // Verify that the transformation isn't quite right.
        assert_abs_diff_eq!(uvfits_xyz, local_xyz, epsilon = 1e0);

        // Now verify cotter's ms XYZ with the constants it uses.
        let cotter_earth_pos = LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        };
        let result = ms_xyz.to_geodetic(cotter_earth_pos);
        assert!(result.is_ok());
        let local_xyz = result.unwrap();
        assert_abs_diff_eq!(uvfits_xyz, local_xyz, epsilon = 1e-6);

        // Now check the conversion of geodetic to geocentric.
        let result = uvfits_xyz.to_geocentric_mwa();
        assert!(result.is_ok());
        let geocentric_xyz = result.unwrap();
        // cotter's MWA coordinates are a little off of what is in mwalib.
        // Verify that the transformation isn't quite right.
        assert_abs_diff_eq!(ms_xyz, geocentric_xyz, epsilon = 1e0);

        // Now verify cotter's ms XYZ with the constants it uses.
        let result = uvfits_xyz.to_geocentric(cotter_earth_pos);
        assert!(result.is_ok());
        let geocentric_xyz = result.unwrap();
        assert_abs_diff_eq!(ms_xyz, geocentric_xyz, epsilon = 1e-6);
    }

    #[test]
    fn xyzs_to_uvws_test() {
        let xyzs = vec![
            XyzGeodetic {
                x: 289.5692922664971,
                y: -585.6749877929688,
                z: -259.3106530519151,
            },
            XyzGeodetic {
                x: 750.5194624923599,
                y: -565.4390258789063,
                z: 665.2348852011041,
            },
        ];
        let phase = HADec::new(6.0163, -0.453121);
        let result: Vec<UVW> = xyzs_to_uvws(&xyzs, phase);
        let expected = UVW {
            u: 102.04605530570603,
            v: -1028.2293398297727,
            w: 0.18220641926160397,
        };
        assert_abs_diff_eq!(result[1], expected, epsilon = 1e-10);
        // Auto-correlations are zero.
        assert_abs_diff_eq!(result[0], UVW::default(),);
        assert_abs_diff_eq!(result[2], UVW::default(),);
    }

    #[test]
    fn xyzs_to_uvws_parallel_test() {
        let xyzs = vec![
            XyzGeodetic {
                x: 289.5692922664971,
                y: -585.6749877929688,
                z: -259.3106530519151,
            },
            XyzGeodetic {
                x: 750.5194624923599,
                y: -565.4390258789063,
                z: 665.2348852011041,
            },
            XyzGeodetic::default(),
            XyzGeodetic::default(),
            XyzGeodetic::default(),
        ];
        let phase = HADec::new(6.0163, -0.453121);
        let serial_result: Vec<UVW> = xyzs_to_uvws(&xyzs, phase);
        let parallel_result: Vec<UVW> = xyzs_to_uvws_parallel(&xyzs, phase);
        assert_eq!(serial_result.len(), 15);
        assert_abs_diff_eq!(Array1::from(serial_result), Array1::from(parallel_result));
    }

    #[test]
    fn xyzs_to_cross_uvws_test() {
        let xyzs = vec![
            XyzGeodetic {
                x: 289.5692922664971,
                y: -585.6749877929688,
                z: -259.3106530519151,
            },
            XyzGeodetic {
                x: 750.5194624923599,
                y: -565.4390258789063,
                z: 665.2348852011041,
            },
        ];
        let phase = HADec::new(6.0163, -0.453121);
        let result: Vec<UVW> = xyzs_to_cross_uvws(&xyzs, phase);
        let expected = UVW {
            u: 102.04605530570603,
            v: -1028.2293398297727,
            w: 0.18220641926160397,
        };
        assert_abs_diff_eq!(
            Array1::from(result),
            Array1::from_elem(1, expected),
            epsilon = 1e-10
        );
    }

    #[test]
    fn xyzs_to_cross_uvws_parallel_test() {
        let xyzs = vec![
            XyzGeodetic {
                x: 289.5692922664971,
                y: -585.6749877929688,
                z: -259.3106530519151,
            },
            XyzGeodetic {
                x: 750.5194624923599,
                y: -565.4390258789063,
                z: 665.2348852011041,
            },
            XyzGeodetic::default(),
            XyzGeodetic::default(),
            XyzGeodetic::default(),
        ];
        let phase = HADec::new(6.0163, -0.453121);
        let serial_result: Vec<UVW> = xyzs_to_cross_uvws(&xyzs, phase);
        let parallel_result: Vec<UVW> = xyzs_to_cross_uvws_parallel(&xyzs, phase);
        assert_eq!(serial_result.len(), 10);
        assert_abs_diff_eq!(Array1::from(serial_result), Array1::from(parallel_result));
    }
}
