// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Conventions of instruments and software packages.
//!
//! Marlu (and `hyperdrive`, `Birli`, `cotter` and the `RTS`) grew up with the
//! `MWA`, whose X dipoles point east-west and whose baselines are antenna 1
//! minus antenna 2. The `IAU`, and the casacore-based software most other
//! instruments use (`CASA`, `DP3`/`LOFAR`, `WSClean`, `ASKAPsoft`), define X
//! as north-south. casacore's UVW calculators (`derivedmscal`, and `DP3`'s)
//! make a baseline antenna 2 minus antenna 1, although the Measurement Sets
//! `CASA` itself writes (`importevla`, `importuvfits`), like cotter's and
//! Birli's, have antenna 1 minus antenna 2, which is also what pyuvdata
//! assumes of any Measurement Set or uvfits file it did not write. casacore's
//! J2000 `UVW`s also differ from Marlu's precession by the equation of the
//! equinoxes and the annual aberration of the phase centre (about 20 arcsec
//! each). These types name the conventions so that code can convert between
//! them and files can record them.
//!
//! The polarisation side is factored the way the
//! [UVH5 format](https://github.com/RadioAstronomySoftwareGroup/pyuvdata/blob/main/docs/references/uvh5_memo.pdf)
//! records it: an [`XOrientation`] (UVH5 `x_orientation`, "east" or "north")
//! and a [`StokesConvention`] (UVH5 `pol_convention`, "avg" for
//! I = (XX+YY)/2 or "sum" for I = XX+YY). A Measurement Set says the same
//! thing natively through the `FEED` table's `RECEPTOR_ANGLE` (east is X at
//! 90°), a uvfits file through the `AIPS AN` table's `POLAA`/`POLAB`, which
//! is also how pyuvdata records (and derives `x_orientation` from) its
//! `feed_angle`. Neither has a slot for "sum", so writers add the keyword
//! pyuvdata uses for it (`pyuvdata_polconv` on an MS main table, `POLCONV`
//! in a uvfits primary header, written only for calibrated [`VisUnits`] as
//! pyuvdata requires), and a keyword for the [`UvwFrame`], which no format
//! records (`marlu_uvw_frame` on an MS main table; `UVWFRAME` in a uvfits
//! primary header, 8 characters so that it is a standard FITS keyword that
//! pyuvdata keeps when it rewrites a file). See
//! [`crate::precession::precess_time_in_frame`] for the geometry behind
//! [`UvwFrame`].

use std::{f64::consts::FRAC_PI_2, fmt::Display, str::FromStr};

use thiserror::Error;

use crate::{Complex, Jones};

/// An unrecognised convention name.
#[derive(Error, Debug, Clone, PartialEq, Eq)]
#[error("'{got}' is not a known {kind}; expected one of {expected}")]
pub struct UnknownConvention {
    kind: &'static str,
    got: String,
    expected: &'static str,
}

/// The orientation of the X dipole/feed (UVH5 `x_orientation`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "lowercase"))]
pub enum XOrientation {
    /// X is east-west (the `MWA`): XX = I - Q, XY = U - iV, YX = U + iV,
    /// YY = I + Q.
    #[default]
    East,
    /// X is north-south (the `IAU` definition, used by casacore/`CASA`, `DP3`,
    /// `WSClean`, `ASKAP`): XX = I + Q, XY = U + iV, YX = U - iV,
    /// YY = I - Q.
    North,
}

impl XOrientation {
    /// The names [`XOrientation::from_str`] accepts (case-insensitive).
    pub const NAMES: &'static str = "east, north";

    /// Position angles (radians, north through east) of the X and Y feeds, as
    /// a Measurement Set's `FEED.RECEPTOR_ANGLE` or a uvfits `AIPS AN` table's
    /// `POLAA`/`POLAB` record them.
    pub fn feed_angles_rad(self) -> [f64; 2] {
        match self {
            XOrientation::East => [FRAC_PI_2, 0.0],
            XOrientation::North => [0.0, FRAC_PI_2],
        }
    }

    /// The orientation recorded by feed angles (see
    /// [`XOrientation::feed_angles_rad`]): X within 45° of east is
    /// [`XOrientation::East`], anything else [`XOrientation::North`]. This is
    /// looser than pyuvdata, which calls feeds "east" or "north" only when
    /// both X and Y angles match exactly, and otherwise has no orientation.
    pub fn from_feed_angle_rad(x_angle_rad: f64) -> XOrientation {
        let a = x_angle_rad.rem_euclid(std::f64::consts::PI);
        if (a - FRAC_PI_2).abs() < std::f64::consts::FRAC_PI_4 {
            XOrientation::East
        } else {
            XOrientation::North
        }
    }
}

impl FromStr for XOrientation {
    type Err = UnknownConvention;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_lowercase().as_str() {
            "east" | "e" => Ok(XOrientation::East),
            "north" | "n" => Ok(XOrientation::North),
            _ => Err(UnknownConvention {
                kind: "x orientation",
                got: s.to_string(),
                expected: Self::NAMES,
            }),
        }
    }
}

impl Display for XOrientation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            XOrientation::East => "east",
            XOrientation::North => "north",
        })
    }
}

/// How instrumental polarisations combine into Stokes parameters (UVH5
/// `pol_convention`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "lowercase"))]
pub enum StokesConvention {
    /// I = (XX + YY) / 2 (the `IAU`, the `MWA`, casacore, ...).
    #[default]
    Avg,
    /// I = XX + YY (`ASKAPsoft`).
    Sum,
}

impl StokesConvention {
    /// The names [`StokesConvention::from_str`] accepts (case-insensitive).
    pub const NAMES: &'static str = "avg, sum";
}

impl FromStr for StokesConvention {
    type Err = UnknownConvention;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_lowercase().as_str() {
            "avg" | "average" | "mean" => Ok(StokesConvention::Avg),
            "sum" => Ok(StokesConvention::Sum),
            _ => Err(UnknownConvention {
                kind: "Stokes convention",
                got: s.to_string(),
                expected: Self::NAMES,
            }),
        }
    }
}

impl Display for StokesConvention {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            StokesConvention::Avg => "avg",
            StokesConvention::Sum => "sum",
        })
    }
}

/// The units of the visibilities (uvfits `BUNIT`, UVH5 `vis_units`).
///
/// pyuvdata reads a Stokes convention (`POLCONV`) only for calibrated data,
/// and refuses a file that records one with uncalibrated visibilities, so
/// writers record [`StokesConvention`] only when the units are calibrated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "lowercase"))]
pub enum VisUnits {
    /// Correlator units, e.g. raw or preprocessed data (uvfits `UNCALIB`).
    #[default]
    Uncalibrated,
    /// Calibrated or modelled flux density (uvfits `JY`).
    Jansky,
}

impl VisUnits {
    /// The names [`VisUnits::from_str`] accepts (case-insensitive).
    pub const NAMES: &'static str = "uncalib, jy";

    /// Whether the visibilities are calibrated.
    pub fn is_calibrated(self) -> bool {
        matches!(self, VisUnits::Jansky)
    }

    /// The uvfits `BUNIT` value (as pyuvdata writes it).
    pub fn bunit(self) -> &'static str {
        match self {
            VisUnits::Uncalibrated => "UNCALIB",
            VisUnits::Jansky => "JY",
        }
    }
}

impl FromStr for VisUnits {
    type Err = UnknownConvention;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_lowercase().as_str() {
            "uncalib" | "uncalibrated" => Ok(VisUnits::Uncalibrated),
            "jy" | "jansky" => Ok(VisUnits::Jansky),
            _ => Err(UnknownConvention {
                kind: "visibility unit",
                got: s.to_string(),
                expected: Self::NAMES,
            }),
        }
    }
}

impl Display for VisUnits {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.bunit())
    }
}

/// How Stokes I, Q, U, V map onto the instrumental polarisations XX, XY, YX,
/// YY: an [`XOrientation`] and a [`StokesConvention`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PolConvention {
    pub x_orientation: XOrientation,
    pub stokes: StokesConvention,
}

impl PolConvention {
    /// The `MWA`: X east-west, I = (XX + YY) / 2.
    pub const MWA: PolConvention = PolConvention {
        x_orientation: XOrientation::East,
        stokes: StokesConvention::Avg,
    };
    /// The `IAU` (casacore/`CASA`, `DP3`/`LOFAR`, `WSClean`): X
    /// north-south, I = (XX + YY) / 2.
    pub const IAU: PolConvention = PolConvention {
        x_orientation: XOrientation::North,
        stokes: StokesConvention::Avg,
    };
    /// `ASKAP` (`ASKAPsoft`): X north-south, I = XX + YY.
    pub const ASKAP: PolConvention = PolConvention {
        x_orientation: XOrientation::North,
        stokes: StokesConvention::Sum,
    };

    /// The names [`PolConvention::from_str`] accepts (case-insensitive):
    /// presets, or `<x orientation>` (avg) or `<x orientation>/<stokes>`.
    pub const NAMES: &'static str = "mwa (alias: oskar), iau (aliases: lofar, casacore, wsclean), askap (alias: askapsoft), or east|north[/avg|sum]";

    /// Convert Stokes flux densities into a [`Jones`] matrix of instrumental
    /// flux densities ordered XX, XY, YX, YY.
    pub fn stokes_to_instrumental(self, i: f64, q: f64, u: f64, v: f64) -> Jones<f64> {
        let scale = match self.stokes {
            StokesConvention::Avg => 1.0,
            StokesConvention::Sum => 0.5,
        };
        match self.x_orientation {
            XOrientation::East => Jones::from([
                Complex::new(scale * (i - q), 0.0),
                Complex::new(scale * u, -scale * v),
                Complex::new(scale * u, scale * v),
                Complex::new(scale * (i + q), 0.0),
            ]),
            XOrientation::North => Jones::from([
                Complex::new(scale * (i + q), 0.0),
                Complex::new(scale * u, scale * v),
                Complex::new(scale * u, -scale * v),
                Complex::new(scale * (i - q), 0.0),
            ]),
        }
    }

    /// Stokes I, Q, U, V from instrumental flux densities ordered XX, XY, YX,
    /// YY (the inverse of [`PolConvention::stokes_to_instrumental`]).
    pub fn instrumental_to_stokes(self, j: Jones<f64>) -> [f64; 4] {
        let (xx, xy, yx, yy) = (j[0], j[1], j[2], j[3]);
        let scale = match self.stokes {
            StokesConvention::Avg => 0.5,
            StokesConvention::Sum => 1.0,
        };
        let (i, q, u, v) = match self.x_orientation {
            XOrientation::East => ((xx + yy).re, (yy - xx).re, (xy + yx).re, (yx - xy).im),
            XOrientation::North => ((xx + yy).re, (xx - yy).re, (xy + yx).re, (xy - yx).im),
        };
        [scale * i, scale * q, scale * u, scale * v]
    }

    /// Position angles (radians) of the X and Y feeds; see
    /// [`XOrientation::feed_angles_rad`].
    pub fn feed_angles_rad(self) -> [f64; 2] {
        self.x_orientation.feed_angles_rad()
    }

    /// A one-line description of the mapping.
    pub fn describe(self) -> String {
        let (preset, mapping) = match self.x_orientation {
            XOrientation::East => ("X east-west", "XX = I-Q, XY = U-iV, YX = U+iV, YY = I+Q"),
            XOrientation::North => ("X north-south", "XX = I+Q, XY = U+iV, YX = U-iV, YY = I-Q"),
        };
        match self.stokes {
            StokesConvention::Avg => format!("{preset}, {mapping}"),
            StokesConvention::Sum => format!("{preset}, half of {mapping} (I = XX+YY)"),
        }
    }
}

impl FromStr for PolConvention {
    type Err = UnknownConvention;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let lower = s.trim().to_lowercase();
        match lower.as_str() {
            // OSKAR's X element lies along its station x (east) axis.
            "mwa" | "oskar" => return Ok(PolConvention::MWA),
            "iau" | "lofar" | "casacore" | "wsclean" => return Ok(PolConvention::IAU),
            "askap" | "askapsoft" => return Ok(PolConvention::ASKAP),
            _ => (),
        }
        let unknown = || UnknownConvention {
            kind: "polarisation convention",
            got: s.to_string(),
            expected: Self::NAMES,
        };
        let (x, stokes) = match lower.split_once('/') {
            Some((x, stokes)) => (x, stokes),
            None => (lower.as_str(), "avg"),
        };
        Ok(PolConvention {
            x_orientation: x.parse().map_err(|_| unknown())?,
            stokes: stokes.parse().map_err(|_| unknown())?,
        })
    }
}

impl Display for PolConvention {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            PolConvention::MWA => f.write_str("mwa"),
            PolConvention::IAU => f.write_str("iau"),
            PolConvention::ASKAP => f.write_str("askap"),
            PolConvention {
                x_orientation,
                stokes,
            } => write!(f, "{x_orientation}/{stokes}"),
        }
    }
}

/// Which J2000 frame an array is placed in, i.e. what the UVWs (and so the
/// visibility phases) mean. See [`crate::precession::precess_time_in_frame`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "lowercase"))]
pub enum UvwFrame {
    /// Marlu's (and `hyperdrive`'s, `Birli`'s and `cotter`'s) frame: the
    /// array is precessed and nutated to J2000 and rotated with the mean
    /// sidereal time, and a baseline is antenna 1 minus antenna 2. This is
    /// also the `AIPS` uvfits convention.
    #[default]
    Hyperdrive,

    /// The frame of the UVWs casacore computes (`derivedmscal`'s
    /// `UVWJ2000`, `measures.to_uvw`) and `DP3` uses, with the baseline sign
    /// of a UVH5 `uvw_array`: apparent sidereal time, the J2000 baselines
    /// rotated to undo the annual aberration of the phase centre, and a
    /// baseline is antenna 2 minus antenna 1 (so visibilities are the complex
    /// conjugates of the hyperdrive frame's). The Measurement Sets `CASA`
    /// itself writes have the opposite baseline sign; Marlu's writers record
    /// this one in `pyuvdata_flip_conj`, and uvfits files are written with
    /// antenna 1 minus antenna 2 baselines (and conjugated visibilities)
    /// whatever the frame.
    #[cfg_attr(
        feature = "serde",
        serde(
            alias = "lofar",
            alias = "casa",
            alias = "wsclean",
            alias = "askap",
            alias = "uvh5",
            alias = "pyuvdata"
        )
    )]
    Casacore,

    /// The frame of an `OSKAR` simulation (with its default
    /// `use_casa_phase_convention`): nothing is precessed - the sky
    /// coordinates are taken as apparent coordinates of date and the array
    /// is rotated with the apparent sidereal time - and a baseline is
    /// antenna 1 minus antenna 2, as in the hyperdrive frame.
    Oskar,
}

impl UvwFrame {
    /// The names [`UvwFrame::from_str`] accepts (case-insensitive).
    pub const NAMES: &'static str = "hyperdrive (aliases: marlu, mwa, aips), casacore (aliases: lofar, casa, wsclean, askap, uvh5, pyuvdata), oskar";

    /// Sign of a baseline vector in this frame relative to Marlu's
    /// (antenna 1 minus antenna 2) convention.
    pub fn baseline_sign(self) -> f64 {
        match self {
            UvwFrame::Hyperdrive | UvwFrame::Oskar => 1.0,
            UvwFrame::Casacore => -1.0,
        }
    }

    /// Whether the array is precessed to J2000 in this frame.
    pub fn precesses(self) -> bool {
        !matches!(self, UvwFrame::Oskar)
    }

    /// A one-line description of the frame.
    pub fn describe(self) -> &'static str {
        match self {
            UvwFrame::Hyperdrive => "hyperdrive (J2000, mean sidereal time, antenna1 - antenna2)",
            UvwFrame::Casacore => {
                "casacore (J2000, apparent sidereal time, phase-centre aberration rotation, antenna2 - antenna1)"
            }
            UvwFrame::Oskar => {
                "oskar (no precession: coordinates of date, apparent sidereal time, antenna1 - antenna2)"
            }
        }
    }
}

impl FromStr for UvwFrame {
    type Err = UnknownConvention;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_lowercase().as_str() {
            "hyperdrive" | "marlu" | "mwa" | "aips" => Ok(UvwFrame::Hyperdrive),
            "casacore" | "lofar" | "casa" | "wsclean" | "askap" | "uvh5" | "pyuvdata" => {
                Ok(UvwFrame::Casacore)
            }
            "oskar" => Ok(UvwFrame::Oskar),
            _ => Err(UnknownConvention {
                kind: "UVW frame",
                got: s.to_string(),
                expected: Self::NAMES,
            }),
        }
    }
}

impl Display for UvwFrame {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            UvwFrame::Hyperdrive => "hyperdrive",
            UvwFrame::Casacore => "casacore",
            UvwFrame::Oskar => "oskar",
        })
    }
}

/// Presets naming an instrument or software package, which fix both a
/// [`PolConvention`] and a [`UvwFrame`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "lowercase"))]
pub enum Convention {
    /// `MWA` polarisations in the hyperdrive frame.
    #[default]
    Mwa,
    /// `IAU` polarisations in the hyperdrive frame.
    Iau,
    /// `IAU` polarisations in the casacore frame: `DP3` and `LOFAR` data.
    Lofar,
    /// `IAU` polarisations in the casacore frame: Measurement Sets made or
    /// read by casacore, `CASA` and `WSClean`, and UVH5 files.
    #[cfg_attr(
        feature = "serde",
        serde(alias = "casa", alias = "wsclean", alias = "uvh5", alias = "pyuvdata")
    )]
    Casacore,
    /// `ASKAP` polarisations (I = XX + YY) in the casacore frame: `ASKAPsoft`.
    #[cfg_attr(feature = "serde", serde(alias = "askapsoft"))]
    Askap,
    /// `OSKAR`: its X element lies along the station x (east) axis and its
    /// Measurement Sets record `RECEPTOR_ANGLE = [x + pi/2, y]`, so its
    /// polarisations are X east-west, and its UVWs are in the
    /// [`UvwFrame::Oskar`] frame.
    Oskar,
}

impl Convention {
    /// The names [`Convention::from_str`] accepts (case-insensitive).
    pub const NAMES: &'static str = "mwa, iau, lofar, casacore (aliases: casa, wsclean, uvh5, pyuvdata), askap (alias: askapsoft), oskar";

    pub fn pol_convention(self) -> PolConvention {
        match self {
            Convention::Mwa | Convention::Oskar => PolConvention::MWA,
            Convention::Iau | Convention::Lofar | Convention::Casacore => PolConvention::IAU,
            Convention::Askap => PolConvention::ASKAP,
        }
    }

    pub fn uvw_frame(self) -> UvwFrame {
        match self {
            Convention::Mwa | Convention::Iau => UvwFrame::Hyperdrive,
            Convention::Lofar | Convention::Casacore | Convention::Askap => UvwFrame::Casacore,
            Convention::Oskar => UvwFrame::Oskar,
        }
    }
}

impl FromStr for Convention {
    type Err = UnknownConvention;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_lowercase().as_str() {
            "mwa" => Ok(Convention::Mwa),
            "iau" => Ok(Convention::Iau),
            "lofar" => Ok(Convention::Lofar),
            "casacore" | "casa" | "wsclean" | "uvh5" | "pyuvdata" => Ok(Convention::Casacore),
            "askap" | "askapsoft" => Ok(Convention::Askap),
            "oskar" => Ok(Convention::Oskar),
            _ => Err(UnknownConvention {
                kind: "convention",
                got: s.to_string(),
                expected: Self::NAMES,
            }),
        }
    }
}

impl Display for Convention {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Convention::Mwa => "mwa",
            Convention::Iau => "iau",
            Convention::Lofar => "lofar",
            Convention::Casacore => "casacore",
            Convention::Askap => "askap",
            Convention::Oskar => "oskar",
        })
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    const I: f64 = 10.0;
    const Q: f64 = 1.0;
    const U: f64 = 2.0;
    const V: f64 = 3.0;

    #[test]
    fn mwa_mapping() {
        let j = PolConvention::MWA.stokes_to_instrumental(I, Q, U, V);
        assert_abs_diff_eq!(j[0], Complex::new(9.0, 0.0));
        assert_abs_diff_eq!(j[1], Complex::new(2.0, -3.0));
        assert_abs_diff_eq!(j[2], Complex::new(2.0, 3.0));
        assert_abs_diff_eq!(j[3], Complex::new(11.0, 0.0));
    }

    #[test]
    fn iau_is_mwa_with_pols_swapped() {
        let m = PolConvention::MWA.stokes_to_instrumental(I, Q, U, V);
        let i = PolConvention::IAU.stokes_to_instrumental(I, Q, U, V);
        assert_abs_diff_eq!(i[0], Complex::new(11.0, 0.0));
        assert_abs_diff_eq!(i[1], Complex::new(2.0, 3.0));
        assert_abs_diff_eq!(i[0], m[3]);
        assert_abs_diff_eq!(i[1], m[2]);
        assert_abs_diff_eq!(i[2], m[1]);
        assert_abs_diff_eq!(i[3], m[0]);
    }

    #[test]
    fn askap_is_half_iau() {
        let a = PolConvention::ASKAP.stokes_to_instrumental(I, Q, U, V);
        let i = PolConvention::IAU.stokes_to_instrumental(I, Q, U, V);
        for p in 0..4 {
            assert_abs_diff_eq!(a[p], i[p] * 0.5);
        }
    }

    #[test]
    fn round_trips() {
        for x in [XOrientation::East, XOrientation::North] {
            for stokes in [StokesConvention::Avg, StokesConvention::Sum] {
                let conv = PolConvention {
                    x_orientation: x,
                    stokes,
                };
                let back = conv.instrumental_to_stokes(conv.stokes_to_instrumental(I, Q, U, V));
                assert_abs_diff_eq!(back[0], I, epsilon = 1e-12);
                assert_abs_diff_eq!(back[1], Q, epsilon = 1e-12);
                assert_abs_diff_eq!(back[2], U, epsilon = 1e-12);
                assert_abs_diff_eq!(back[3], V, epsilon = 1e-12);
                assert_eq!(conv.to_string().parse::<PolConvention>().unwrap(), conv);
            }
        }
    }

    #[test]
    fn feed_angles() {
        let [a, b] = PolConvention::MWA.feed_angles_rad();
        assert_abs_diff_eq!(a, FRAC_PI_2);
        assert_abs_diff_eq!(b, 0.0);
        let [a, b] = PolConvention::IAU.feed_angles_rad();
        assert_abs_diff_eq!(a, 0.0);
        assert_abs_diff_eq!(b, FRAC_PI_2);
        assert_eq!(
            XOrientation::from_feed_angle_rad(FRAC_PI_2),
            XOrientation::East
        );
        assert_eq!(
            XOrientation::from_feed_angle_rad(-FRAC_PI_2),
            XOrientation::East
        );
        assert_eq!(XOrientation::from_feed_angle_rad(0.0), XOrientation::North);
        // LOFAR dipoles are at +-45 degrees but their data are IAU.
        assert_eq!(
            XOrientation::from_feed_angle_rad(-std::f64::consts::FRAC_PI_4),
            XOrientation::North
        );
    }

    #[test]
    fn names_and_aliases() {
        for (name, expected) in [
            ("mwa", PolConvention::MWA),
            ("MWA", PolConvention::MWA),
            ("iau", PolConvention::IAU),
            ("lofar", PolConvention::IAU),
            ("casacore", PolConvention::IAU),
            ("wsclean", PolConvention::IAU),
            ("OSKAR", PolConvention::MWA),
            ("askap", PolConvention::ASKAP),
            ("askapsoft", PolConvention::ASKAP),
            ("east", PolConvention::MWA),
            ("north/avg", PolConvention::IAU),
            ("north/sum", PolConvention::ASKAP),
        ] {
            assert_eq!(name.parse::<PolConvention>().unwrap(), expected, "{name}");
        }
        assert!("circular".parse::<PolConvention>().is_err());
        assert!("east/circular".parse::<PolConvention>().is_err());
        assert_eq!(PolConvention::MWA.to_string(), "mwa");
        assert_eq!(
            PolConvention {
                x_orientation: XOrientation::East,
                stokes: StokesConvention::Sum
            }
            .to_string(),
            "east/sum"
        );
        assert_eq!("lofar".parse::<UvwFrame>().unwrap(), UvwFrame::Casacore);
        assert_eq!("uvh5".parse::<UvwFrame>().unwrap(), UvwFrame::Casacore);
        assert_eq!("marlu".parse::<UvwFrame>().unwrap(), UvwFrame::Hyperdrive);
        assert_eq!(
            "wsclean".parse::<Convention>().unwrap(),
            Convention::Casacore
        );
        for c in [
            Convention::Mwa,
            Convention::Iau,
            Convention::Lofar,
            Convention::Casacore,
            Convention::Askap,
            Convention::Oskar,
        ] {
            assert_eq!(c.to_string().parse::<Convention>().unwrap(), c);
        }
        for f in [UvwFrame::Hyperdrive, UvwFrame::Casacore, UvwFrame::Oskar] {
            assert_eq!(f.to_string().parse::<UvwFrame>().unwrap(), f);
        }
    }

    #[test]
    fn presets() {
        assert_eq!(Convention::Lofar.pol_convention(), PolConvention::IAU);
        assert_eq!(Convention::Lofar.uvw_frame(), UvwFrame::Casacore);
        assert_eq!(Convention::Askap.pol_convention(), PolConvention::ASKAP);
        assert_eq!(Convention::Askap.uvw_frame(), UvwFrame::Casacore);
        assert_eq!(Convention::Oskar.pol_convention(), PolConvention::MWA);
        assert_eq!(Convention::Oskar.uvw_frame(), UvwFrame::Oskar);
        assert_eq!("oskar".parse::<UvwFrame>().unwrap(), UvwFrame::Oskar);
        assert!(!UvwFrame::Oskar.precesses() && UvwFrame::Casacore.precesses());
        assert_eq!(Convention::Mwa.uvw_frame(), UvwFrame::Hyperdrive);
        assert_abs_diff_eq!(UvwFrame::Casacore.baseline_sign(), -1.0);
        assert_abs_diff_eq!(UvwFrame::Hyperdrive.baseline_sign(), 1.0);
        assert_abs_diff_eq!(UvwFrame::Oskar.baseline_sign(), 1.0);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn serde_names_and_aliases() {
        let x: XOrientation = serde_json::from_str("\"north\"").unwrap();
        assert_eq!(x, XOrientation::North);
        let s: StokesConvention = serde_json::from_str("\"sum\"").unwrap();
        assert_eq!(s, StokesConvention::Sum);
        let p: PolConvention =
            serde_json::from_str("{\"x_orientation\": \"east\", \"stokes\": \"avg\"}").unwrap();
        assert_eq!(p, PolConvention::MWA);
        let f: UvwFrame = serde_json::from_str("\"askap\"").unwrap();
        assert_eq!(f, UvwFrame::Casacore);
        let c: Convention = serde_json::from_str("\"casa\"").unwrap();
        assert_eq!(c, Convention::Casacore);
        assert_eq!(
            serde_json::to_string(&UvwFrame::Casacore).unwrap(),
            "\"casacore\""
        );
    }
}
