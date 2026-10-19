// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Code to handle precession (including nutation).
//!
//! This is a really good resource:
//! <https://lweb.cfa.harvard.edu/~jzhao/times.html>
//!
//! Some easy-to-read info is here:
//! <https://en.wikipedia.org/wiki/Astronomical_nutation>
//!
//! A harder-to-read source of info is here:
//! <https://www.aanda.org/articles/aa/pdf/2003/48/aa4068.pdf>

use std::f64::consts::TAU;

use erfa::constants::{ERFA_CMPS, ERFA_DAU, ERFA_DAYSEC, ERFA_DJM0};
use hifitime::{Duration, Epoch};

use crate::{pal, HADec, RADec, XyzGeodetic};

#[derive(Debug, Clone, Copy)]
pub struct PrecessionInfo {
    /// Bias procession rotation matrix.
    rotation_matrix: [[f64; 3]; 3],

    /// The phase centre's hour angle (against [`PrecessionInfo::lmst_j2000`])
    /// and declination in the J2000 frame of [`PrecessionInfo::precess_xyz`],
    /// i.e. `phase_centre.to_hadec(lmst_j2000)`. Project precessed positions
    /// onto this (or any other J2000 direction with `to_hadec(lmst_j2000)`)
    /// to get `UVWs`; the annual aberration is already in the positions.
    pub hadec_j2000: HADec,

    /// The LMST of the current epoch.
    pub lmst: f64,

    /// The precessed LMST in the J2000 epoch.
    pub lmst_j2000: f64,

    /// The precessed array latitude in the J2000 epoch.
    pub array_latitude_j2000: f64,
}

impl PrecessionInfo {
    // Blatently stolen from cotter.
    pub fn precess_xyz(&self, xyzs: &[XyzGeodetic]) -> Vec<XyzGeodetic> {
        let (sep, cep) = self.lmst.sin_cos();
        let (s2000, c2000) = self.lmst_j2000.sin_cos();

        xyzs.iter()
            .map(|xyz| {
                // rotate to frame with x axis at zero RA
                let xpr = cep * xyz.x - sep * xyz.y;
                let ypr = sep * xyz.x + cep * xyz.y;
                let zpr = xyz.z;

                let rmat = &self.rotation_matrix;
                let xpr2 = (rmat[0][0]) * xpr + (rmat[0][1]) * ypr + (rmat[0][2]) * zpr;
                let ypr2 = (rmat[1][0]) * xpr + (rmat[1][1]) * ypr + (rmat[1][2]) * zpr;
                let zpr2 = (rmat[2][0]) * xpr + (rmat[2][1]) * ypr + (rmat[2][2]) * zpr;

                // rotate back to frame with xp pointing out at lmst2000
                XyzGeodetic {
                    x: c2000 * xpr2 + s2000 * ypr2,
                    y: -s2000 * xpr2 + c2000 * ypr2,
                    z: zpr2,
                }
            })
            .collect()
    }

    #[deprecated = "use `PrecessionInfo::precess_xyz` instead"]
    pub fn precess_xyz_parallel(&self, xyzs: &[XyzGeodetic]) -> Vec<XyzGeodetic> {
        self.precess_xyz(xyzs)
    }
}

/// Get the local mean sidereal time. `time` should be in the UTC frame, and
/// `dut1` (i.e. UT1 - UTC) provides a better estimate of the LMST. If DUT1
/// isn't known, then a [`Duration`] of 0 seconds can be used; the results are
/// wrong by up to 0.9 seconds.
pub fn get_lmst(array_longitude_rad: f64, time: Epoch, dut1: Duration) -> f64 {
    let ut1 = (time + dut1).to_mjd_utc_days();
    let gmst = pal::palGmst(ut1);
    (gmst + array_longitude_rad) % TAU
}

/// Obtain precessed coordinate information. `time` should be in the UTC frame,
/// and `dut1` (i.e. UT1 - UTC) provides a better estimate of the LMST. If DUT1
/// isn't known, then a [`Duration`] of 0 seconds can be used; the results are
/// wrong by up to 0.9 seconds.
///
/// This function is very similar to `cotter`'s `PrepareTimestepUVW`, with one
/// difference. `cotter` (and earlier versions of Marlu) projected the precessed
/// positions onto the *aberrated* phase centre's hour angle and declination,
/// which gets `w` right but rolls `u` and `v` by a few arcsec (the RA shift
/// times the sine of the declination), and leaves the J2000 positions
/// inconsistent with `phase_centre.to_hadec(lmst_j2000)`, which is how
/// modellers (e.g. `hyperdrive`) project them. Here, as in `casacore`'s
/// `MBaseline` conversion and `pyuvdata`, the annual aberration of the phase
/// centre is instead taken out of the positions by the minimal rigid rotation
/// (`casacore`'s `MeasMath::deapplyAberration`), so that `UVWs` are
/// projected onto the catalogue phase centre: `w` is unchanged, and a model
/// built from catalogue source positions and these `UVWs` sees the sky where
/// the (aberrated) data have it. `cotter` also rotated the array by the mean
/// sidereal time; the array is rotated here by the apparent one, which
/// matches the nutation in the precession matrix.
pub fn precess_time(
    array_longitude_rad: f64,
    array_latitude_rad: f64,
    phase_centre: RADec,
    time: Epoch,
    dut1: Duration,
) -> PrecessionInfo {
    let lmst = get_lmst(array_longitude_rad, time, dut1);
    // The matrix below includes nutation, so it takes positions from the true
    // equator and equinox of date: the array has to be rotated by the
    // apparent sidereal time, not the mean one. The difference is the
    // equation of the equinoxes (up to ~17 arcsec).
    let ut1 = (time + dut1).to_mjd_utc_days();
    let tt = time.to_mjd_tt_days();
    let eqeq = erfa::time::gst06a(ERFA_DJM0, ut1, ERFA_DJM0, tt)
        - erfa::time::gmst06(ERFA_DJM0, ut1, ERFA_DJM0, tt);
    let last = lmst + eqeq;

    let j2000 = 2000.0;
    let mjd = (time + dut1).to_mjd_utc_days();
    let mut rotation_matrix = [[0.0; 3]; 3];
    unsafe { pal::palPrenut(j2000, mjd, rotation_matrix.as_mut_ptr()) };

    // Transpose the rotation matrix.
    let mut rotation_matrix = {
        let mut new = [[0.0; 3]; 3];
        let old = rotation_matrix;
        for (i, old) in old.iter().enumerate() {
            for (j, new) in new.iter_mut().enumerate() {
                new[i] = old[j];
            }
        }
        new
    };

    // The zenith gives the J2000 sidereal time and latitude; these don't
    // include aberration.
    let (lmst_j2000, array_latitude_j2000) =
        rotate_radec(&mut rotation_matrix, last, array_latitude_rad);

    // `precess_xyz` rotates positions by `lmst`; the rest of the way to the
    // apparent sidereal time goes in the matrix.
    let (s, c) = eqeq.sin_cos();
    let eqeq_rotation = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
    let rotation_matrix = erfa::vectors_and_matrices::multiply_matrices(
        deaberration(phase_centre, time),
        erfa::vectors_and_matrices::multiply_matrices(rotation_matrix, eqeq_rotation),
    );

    PrecessionInfo {
        rotation_matrix,
        hadec_j2000: phase_centre.to_hadec(lmst_j2000),
        lmst,
        lmst_j2000,
        array_latitude_j2000,
    }
}

/// The rigid rotation taking the annual aberration of `phase_centre` out of
/// J2000 positions: the minimal rotation from the catalogue phase centre to
/// the direction whose aberrated position is the catalogue one (`casacore`'s
/// `MeasMath::deapplyAberration` with `rotateShift`).
fn deaberration(phase_centre: RADec, time: Epoch) -> [[f64; 3]; 3] {
    let tt = time.to_mjd_tt_days();
    let (_, _, pvb) = erfa::earth::position_velocity_00(ERFA_DJM0, tt);
    let v_over_c = pvb[1].map(|v| v * ERFA_DAU / ERFA_DAYSEC / ERFA_CMPS);
    let p_cat = erfa::transform::spherical_to_cartesian(phase_centre.ra, phase_centre.dec);
    let p_geo = deaberrate(p_cat, v_over_c);
    rotation_taking(p_cat, p_geo)
}

/// Relativistic annual aberration of unit direction `p` (`casacore`
/// `MeasMath::applyAberration`).
fn aberrate(p: [f64; 3], v: [f64; 3]) -> [f64; 3] {
    let g1 = (1.0 - dot(v, v)).sqrt();
    let g2 = dot(p, v);
    let a = (g1 - 1.0 - g2) / (1.0 + g2);
    let b = (1.0 + g2 / (1.0 + g1)) / (1.0 + g2);
    normalise([
        p[0] + a * p[0] + b * v[0],
        p[1] + a * p[1] + b * v[1],
        p[2] + a * p[2] + b * v[2],
    ])
}

/// The direction whose aberrated position is `p_app` (`casacore`
/// `MeasMath::deapplyAberration`).
fn deaberrate(p_app: [f64; 3], v: [f64; 3]) -> [f64; 3] {
    let mut p = normalise([p_app[0] - v[0], p_app[1] - v[1], p_app[2] - v[2]]);
    for _ in 0..50 {
        let q = aberrate(p, v);
        let delta = [q[0] - p_app[0], q[1] - p_app[1], q[2] - p_app[2]];
        if delta.iter().all(|d| d.abs() < 1e-14) {
            break;
        }
        p = normalise([p[0] - delta[0], p[1] - delta[1], p[2] - delta[2]]);
    }
    p
}

/// Minimal rigid rotation taking unit vector `a` to unit vector `b`.
fn rotation_taking(a: [f64; 3], b: [f64; 3]) -> [[f64; 3]; 3] {
    let axis = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let s = dot(axis, axis).sqrt();
    let c = dot(a, b);
    let mut r = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    if s < 1e-18 {
        return r;
    }
    let k = [axis[0] / s, axis[1] / s, axis[2] / s];
    let kx = [[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]];
    let kx2 = erfa::vectors_and_matrices::multiply_matrices(kx, kx);
    for (i, row) in r.iter_mut().enumerate() {
        for (j, e) in row.iter_mut().enumerate() {
            *e += s * kx[i][j] + (1.0 - c) * kx2[i][j];
        }
    }
    r
}

#[cfg(test)]
fn transpose(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0; 3]; 3];
    for (i, row) in m.iter().enumerate() {
        for (j, e) in row.iter().enumerate() {
            t[j][i] = *e;
        }
    }
    t
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn normalise(v: [f64; 3]) -> [f64; 3] {
    let n = dot(v, v).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}

// Blatently stolen from cotter.
fn rotate_radec(rotation_matrix: &mut [[f64; 3]; 3], ra: f64, dec: f64) -> (f64, f64) {
    let mut v1 = [0.0; 3];
    let mut v2 = [0.0; 3];

    unsafe {
        pal::palDcs2c(ra, dec, v1.as_mut_ptr());
        pal::palDmxv(
            rotation_matrix.as_mut_ptr(),
            v1.as_mut_ptr(),
            v2.as_mut_ptr(),
        );
        let mut ra2 = 0.0;
        let mut dec2 = 0.0;
        pal::palDcc2s(v2.as_mut_ptr(), &mut ra2, &mut dec2);
        ra2 = pal::palDranrm(ra2);
        (ra2, dec2)
    }
}

#[cfg(test)]
mod tests {
    use approx::{assert_abs_diff_eq, assert_abs_diff_ne};
    use hifitime::Unit;
    use std::str::FromStr;

    use super::*;
    use crate::constants::{MWA_LAT_RAD, MWA_LONG_RAD};

    // Expected values are taken from astropy 5.0.4, calculated with e.g.
    //
    // loc = EarthLocation(lat=-0.4660608448386394*u.rad, lon=2.0362898668561042*u.rad, height=377.827*u.m)
    // np.deg2rad(Time("1090008642", format="gps", scale="utc", location=loc).sidereal_time("mean").value*15)
    //
    // DUT1 values are calculated with e.g.
    //
    // Time("1090008643", format="gps", scale="utc", location=loc).delta_ut1_utc
    #[test]
    fn aberration_round_trip() {
        let v = [1e-4, -5e-5, 2e-5];
        let p = normalise([0.3, -0.5, 0.8]);
        let q = aberrate(p, v);
        let back = deaberrate(q, v);
        for i in 0..3 {
            assert_abs_diff_eq!(back[i], p[i], epsilon = 1e-13);
        }
        let moved = erfa::vectors_and_matrices::mat_mul_pvec(rotation_taking(p, q), p);
        for i in 0..3 {
            assert_abs_diff_eq!(moved[i], q[i], epsilon = 1e-13);
        }
    }

    #[test]
    fn uvws_are_projected_onto_the_catalogue_phase_centre() {
        use crate::pos::xyz::xyzs_to_cross_uvws;

        // 1254670392 (EoR0, October 2019).
        let phase_centre = RADec::from_degrees(0.0, -27.0);
        let time = Epoch::from_gpst_seconds(1254670392.0);
        let p = precess_time(
            MWA_LONG_RAD,
            MWA_LAT_RAD,
            phase_centre,
            time,
            Duration::default(),
        );
        // What modellers project onto.
        let pc_hadec = phase_centre.to_hadec(p.lmst_j2000);
        assert_abs_diff_eq!(p.hadec_j2000.ha, pc_hadec.ha);
        assert_abs_diff_eq!(p.hadec_j2000.dec, pc_hadec.dec);

        // cotter's frame: positions without the de-aberration rotation,
        // projected onto the aberrated phase centre.
        let cotter = PrecessionInfo {
            rotation_matrix: erfa::vectors_and_matrices::multiply_matrices(
                transpose(deaberration(phase_centre, time)),
                p.rotation_matrix,
            ),
            ..p
        };
        let (_, _, pvb) = erfa::earth::position_velocity_00(ERFA_DJM0, time.to_mjd_tt_days());
        let v = pvb[1].map(|v| v * ERFA_DAU / ERFA_DAYSEC / ERFA_CMPS);
        let (ra_aber, dec_aber) = erfa::transform::cartesian_to_spherical(aberrate(
            erfa::transform::spherical_to_cartesian(phase_centre.ra, phase_centre.dec),
            v,
        ));
        let cotter_hadec = HADec::from_radians(p.lmst_j2000 - ra_aber, dec_aber);

        let xyzs = [
            XyzGeodetic::default(),
            XyzGeodetic {
                x: 3000.0,
                y: 0.0,
                z: 0.0,
            },
            XyzGeodetic {
                x: 0.0,
                y: 3000.0,
                z: 0.0,
            },
            XyzGeodetic {
                x: 0.0,
                y: 0.0,
                z: 3000.0,
            },
            XyzGeodetic {
                x: 1500.0,
                y: -2000.0,
                z: 800.0,
            },
        ];
        let new = xyzs_to_cross_uvws(&p.precess_xyz(&xyzs), p.hadec_j2000);
        let old = xyzs_to_cross_uvws(&cotter.precess_xyz(&xyzs), cotter_hadec);

        // The aberration is the same rigid rotation in both, so `w` and the
        // baseline lengths agree; `u` and `v` differ by a roll about the phase
        // centre of the RA aberration times sin(dec): the minimal rotation
        // doesn't keep north pointing at the pole, the RA/Dec shift does.
        let d_ra = (ra_aber - phase_centre.ra + std::f64::consts::PI).rem_euclid(TAU)
            - std::f64::consts::PI;
        let expected_roll = d_ra * phase_centre.dec.sin();
        assert!(expected_roll.abs().to_degrees() * 3600.0 > 5.0);
        for (n, o) in new.iter().zip(old.iter()) {
            assert_abs_diff_eq!(n.w, o.w, epsilon = 1e-4);
            assert_abs_diff_eq!(n.u.hypot(n.v), o.u.hypot(o.v), epsilon = 1e-4);
            let roll = (o.u * n.v - o.v * n.u).atan2(o.u * n.u + o.v * n.v);
            assert_abs_diff_eq!(roll, expected_roll, epsilon = 0.02 * expected_roll.abs());
        }
    }

    // Four stations around the MWA (ITRF, metres; 3-32 km baselines), two
    // timestamps, and the UVWs CALC 11 (the VLBI delay model used by DiFX, via
    // pycalc11, `noatmo` mode, no polar motion, UT1-UTC = -0.151522 s) gives
    // for them towards RA 0, Dec -27 (J2000). CALC's u and v are the gradient
    // of the full relativistic delay over catalogue J2000 offsets, and w is
    // the delay times c, so they don't depend on a projection convention.
    const CALC_ITRF: [[f64; 3]; 4] = [
        [-2559454.0790495784, 5095372.143311582, -2849057.185786278],
        [-2562134.8794759056, 5094025.551660707, -2849057.185786278],
        [-2560059.1979564168, 5096576.81662318, -2846377.149719187],
        [-2588279.1463356433, 5085921.804508169, -2840123.732229308],
    ];
    const CALC_GPS_TIMES: [f64; 2] = [1254670338.0, 1254670578.0];
    const CALC_UVWS: [(usize, usize, usize, [f64; 3]); 12] = [
        (
            0,
            0,
            1,
            [2975.2923784119776, -174.31783798732067, -343.5356172900647],
        ),
        (
            0,
            0,
            2,
            [173.1940623857081, 2995.060250000395, -19.76255600526929],
        ),
        (
            0,
            0,
            3,
            [30330.23786973243, 8240.354734546112, -3501.2349972864613],
        ),
        (
            0,
            1,
            2,
            [-2802.0983160262695, 3169.378087987716, 323.7730612847954],
        ),
        (
            0,
            1,
            3,
            [27354.945491320454, 8414.672572533433, -3157.6993799963966],
        ),
        (
            0,
            2,
            3,
            [30157.043807346723, 5245.294484545717, -3481.472441281192],
        ),
        (
            1,
            0,
            1,
            [2968.095101361163, -197.84090606482278, -389.91919966600835],
        ),
        (
            1,
            0,
            2,
            [196.56548183760606, 2993.5968511652736, -22.64675517939031],
        ),
        (
            1,
            0,
            3,
            [30336.169577636523, 8000.246257009172, -3974.684804080054],
        ),
        (
            1,
            1,
            2,
            [-2771.529619523557, 3191.4377572300964, 367.27244448661804],
        ),
        (
            1,
            1,
            3,
            [27368.07447627536, 8198.087163073995, -3584.765604414046],
        ),
        (
            1,
            2,
            3,
            [30139.604095798917, 5006.649405843898, -3952.038048900664],
        ),
    ];

    #[test]
    fn uvws_match_calc() {
        use crate::{pos::xyz::xyzs_to_cross_uvws, XyzGeocentric};

        let centre = XyzGeocentric {
            x: CALC_ITRF.iter().map(|p| p[0]).sum::<f64>() / 4.0,
            y: CALC_ITRF.iter().map(|p| p[1]).sum::<f64>() / 4.0,
            z: CALC_ITRF.iter().map(|p| p[2]).sum::<f64>() / 4.0,
        };
        let array_pos = centre.to_earth_wgs84();
        let xyzs: Vec<XyzGeodetic> = CALC_ITRF
            .iter()
            .map(|p| {
                XyzGeocentric {
                    x: p[0],
                    y: p[1],
                    z: p[2],
                }
                .to_geodetic(array_pos)
            })
            .collect();
        let phase_centre = RADec::from_degrees(0.0, -27.0);
        for (i_time, gps) in CALC_GPS_TIMES.iter().enumerate() {
            let time = Epoch::from_gpst_seconds(*gps);
            let p = precess_time(
                array_pos.longitude_rad,
                array_pos.latitude_rad,
                phase_centre,
                time,
                Duration::from_seconds(-0.151522),
            );
            let uvws = xyzs_to_cross_uvws(&p.precess_xyz(&xyzs), p.hadec_j2000);
            // CALC differentiates over catalogue offsets, so its u and v also
            // carry the aberration's compression of the field, a scale of
            // 1 - v.s/c (43 ppm here), which a rigid rotation leaves out.
            let (_, _, pvb) = erfa::earth::position_velocity_00(ERFA_DJM0, time.to_mjd_tt_days());
            let v = pvb[1].map(|v| v * ERFA_DAU / ERFA_DAYSEC / ERFA_CMPS);
            let s = erfa::transform::spherical_to_cartesian(phase_centre.ra, phase_centre.dec);
            let scale = 1.0 - dot(v, s);
            let mut i_bl = 0;
            for a1 in 0..4 {
                for a2 in a1 + 1..4 {
                    let expected = CALC_UVWS
                        .iter()
                        .find(|(t, x, y, _)| *t == i_time && *x == a1 && *y == a2)
                        .unwrap()
                        .3;
                    // CALC's baselines have the opposite sign.
                    let got = uvws[i_bl] * -1.0;
                    let km = (expected[0].powi(2) + expected[1].powi(2) + expected[2].powi(2))
                        .sqrt()
                        / 1e3;
                    // Before taking the aberration out of the positions, u
                    // and v were rolled by 9 arcsec (44 mm per km); with the
                    // mean sidereal time, w was off by 70 mm per km.
                    let duv = (scale * got.u - expected[0]).hypot(scale * got.v - expected[1]);
                    assert!(duv < 2e-3 * km, "baseline {a1}-{a2}: uv off by {duv} m");
                    let dw = (got.w - expected[2]).abs();
                    assert!(dw < 4e-3 * km, "baseline {a1}-{a2}: w off by {dw} m");
                    i_bl += 1;
                }
            }
        }
    }

    #[test]
    fn test_get_lst() {
        let epoch = Epoch::from_gpst_seconds(1090008642.0);
        let dut1 = Duration::from_f64(-0.31295757, Unit::Second);
        assert_abs_diff_eq!(
            get_lmst(MWA_LONG_RAD, epoch, dut1),
            6.262065126600022,
            epsilon = 1e-9
        );

        let epoch = Epoch::from_gpst_seconds(1090008643.0);
        let dut1 = Duration::from_f64(-0.31295757, Unit::Second);
        assert_abs_diff_eq!(
            get_lmst(MWA_LONG_RAD, epoch, dut1),
            6.26213804775838,
            epsilon = 1e-9
        );

        let epoch = Epoch::from_gpst_seconds(1090008647.0);
        let dut1 = Duration::from_f64(-0.31295758, Unit::Second);
        assert_abs_diff_eq!(
            get_lmst(MWA_LONG_RAD, epoch, dut1),
            6.262429732391841,
            epsilon = 1e-9
        );

        let epoch = Epoch::from_gpst_seconds(1090008644.0);
        let dut1 = Duration::from_f64(-0.31295757, Unit::Second);
        assert_abs_diff_eq!(
            get_lmst(MWA_LONG_RAD, epoch, dut1),
            6.262210968916753,
            epsilon = 1e-9
        );
    }

    /// [`PrecessionInfo`]'s matrix without the de-aberration rotation.
    fn precession_matrix(p: &PrecessionInfo, phase_centre: RADec, time: Epoch) -> [[f64; 3]; 3] {
        erfa::vectors_and_matrices::multiply_matrices(
            transpose(deaberration(phase_centre, time)),
            p.rotation_matrix,
        )
    }

    #[test]
    // TODO: reduce cognitive complexity
    #[allow(clippy::cognitive_complexity)]
    fn test_no_precession_at_j2000() {
        // Jack is (pretty) confident that the relatively large errors here are
        // due to nutation, which gets corrected at the same time as precession.

        // https://en.wikipedia.org/wiki/Epoch_(astronomy)#Julian_dates_and_J2000
        let j2000_epoch = Epoch::from_str("2000-01-01T11:58:55.816 UTC").unwrap();

        let phase_centre = RADec::from_degrees(0.0, -27.0);
        let eor0 = precess_time(
            MWA_LONG_RAD,
            MWA_LAT_RAD,
            phase_centre,
            j2000_epoch,
            Duration::from_total_nanoseconds(0),
        );
        let rm0 = precession_matrix(&eor0, phase_centre, j2000_epoch);
        assert_abs_diff_eq!(rm0[0][0], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(rm0[1][1], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(rm0[2][2], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(rm0[0][1], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(rm0[0][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(rm0[1][0], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(rm0[1][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(rm0[2][0], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(rm0[2][1], 0.0, epsilon = 3e-5);

        assert_abs_diff_eq!(eor0.lmst_j2000 - eor0.hadec_j2000.ha, 0.0, epsilon = 1e-4);
        assert_abs_diff_eq!(eor0.array_latitude_j2000, MWA_LAT_RAD, epsilon = 1e-4);
        assert_abs_diff_eq!(eor0.lmst, eor0.lmst_j2000, epsilon = 1e-4);

        assert_abs_diff_eq!(eor0.lmst, 0.6433259676052971, epsilon = 1e-4);

        let phase_centre = RADec::from_degrees(60.0, -30.0);
        let eor1 = precess_time(
            MWA_LONG_RAD,
            MWA_LAT_RAD,
            phase_centre,
            j2000_epoch,
            Duration::from_total_nanoseconds(0),
        );
        let rm1 = precession_matrix(&eor1, phase_centre, j2000_epoch);
        assert_abs_diff_eq!(rm1[0][0], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(rm1[1][1], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(rm1[2][2], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(rm1[0][1], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(rm1[0][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(rm1[1][0], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(rm1[1][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(rm1[2][0], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(rm1[2][1], 0.0, epsilon = 3e-5);

        assert_abs_diff_eq!(
            eor1.lmst_j2000 - eor1.hadec_j2000.ha,
            1.0471975511965976,
            epsilon = 1e-4
        );
        assert_abs_diff_eq!(eor1.array_latitude_j2000, MWA_LAT_RAD, epsilon = 1e-4);
        assert_abs_diff_eq!(eor1.lmst, eor1.lmst_j2000, epsilon = 1e-4);

        assert_abs_diff_eq!(eor1.lmst, 0.6433259676052971, epsilon = 1e-4);
    }

    #[test]
    // TODO: reduce cognitive complexity
    #[allow(clippy::cognitive_complexity)]
    fn test_no_precession_at_j2000_with_dut1() {
        // Note that the test values are exactly the same as above, despite
        // using a DUT1.

        // https://en.wikipedia.org/wiki/Epoch_(astronomy)#Julian_dates_and_J2000
        let j2000_epoch = Epoch::from_str("2000-01-01T11:58:55.816 UTC").unwrap();

        let phase_centre = RADec::from_degrees(0.0, -27.0);
        let eor0 = precess_time(
            MWA_LONG_RAD,
            MWA_LAT_RAD,
            phase_centre,
            j2000_epoch,
            Duration::from_f64(0.35494706, Unit::Second),
        );
        let rm0 = precession_matrix(&eor0, phase_centre, j2000_epoch);
        assert_abs_diff_eq!(rm0[0][0], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(rm0[1][1], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(rm0[2][2], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(rm0[0][1], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(rm0[0][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(rm0[1][0], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(rm0[1][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(rm0[2][0], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(rm0[2][1], 0.0, epsilon = 3e-5);

        assert_abs_diff_eq!(eor0.lmst_j2000 - eor0.hadec_j2000.ha, 0.0, epsilon = 1e-4);
        assert_abs_diff_eq!(eor0.array_latitude_j2000, MWA_LAT_RAD, epsilon = 1e-4);
        assert_abs_diff_eq!(eor0.lmst, eor0.lmst_j2000, epsilon = 1e-4);

        assert_abs_diff_eq!(eor0.lmst, 0.6433259676052971, epsilon = 1e-4);

        let phase_centre = RADec::from_degrees(60.0, -30.0);
        let eor1 = precess_time(
            MWA_LONG_RAD,
            MWA_LAT_RAD,
            phase_centre,
            j2000_epoch,
            Duration::from_f64(0.35494706, Unit::Second),
        );
        let rm1 = precession_matrix(&eor1, phase_centre, j2000_epoch);
        assert_abs_diff_eq!(rm1[0][0], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(rm1[1][1], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(rm1[2][2], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(rm1[0][1], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(rm1[0][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(rm1[1][0], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(rm1[1][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(rm1[2][0], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(rm1[2][1], 0.0, epsilon = 3e-5);

        assert_abs_diff_eq!(
            eor1.lmst_j2000 - eor1.hadec_j2000.ha,
            1.0471975511965976,
            epsilon = 1e-4
        );
        assert_abs_diff_eq!(eor1.array_latitude_j2000, MWA_LAT_RAD, epsilon = 1e-4);
        assert_abs_diff_eq!(eor1.lmst, eor1.lmst_j2000, epsilon = 1e-4);

        assert_abs_diff_eq!(eor1.lmst, 0.6433259676052971, epsilon = 1e-4);
    }

    #[test]
    fn test_precess_1065880128_to_j2000() {
        let epoch = Epoch::from_gpst_seconds(1065880128.0);
        let phase_centre = RADec::from_degrees(0.0, -27.0);

        let p = precess_time(
            MWA_LONG_RAD,
            MWA_LAT_RAD,
            phase_centre,
            epoch,
            Duration::from_total_nanoseconds(0),
        );
        let rm = precession_matrix(&p, phase_centre, epoch);
        // How do I know this is right? Good question! ... I don't.
        assert_abs_diff_eq!(rm[0][0], 1.0, epsilon = 1e-5);
        assert_abs_diff_eq!(rm[1][1], 1.0, epsilon = 1e-5);
        assert_abs_diff_eq!(rm[2][2], 1.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[0][1], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[0][2], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[1][0], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[1][2], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[2][0], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[2][1], 0.0, epsilon = 1e-5);

        // Expected J2000 values from ERFA: the zenith at the apparent sidereal
        // time (`eraGst06a`), taken to J2000 with `eraPnm06a`.
        assert_abs_diff_eq!(p.hadec_j2000.ha, 6.071569090772114, epsilon = 1e-8);
        assert_abs_diff_eq!(p.hadec_j2000.dec, -0.47123889803846897, epsilon = 1e-10);
        assert_abs_diff_eq!(p.lmst, 6.0747789094260245, epsilon = 1e-10);
        assert_abs_diff_eq!(p.lmst_j2000, 6.071569090772114, epsilon = 1e-8);
        assert_abs_diff_eq!(p.array_latitude_j2000, -0.46739656091984, epsilon = 1e-8);
        assert_abs_diff_ne!(p.array_latitude_j2000, MWA_LAT_RAD, epsilon = 1e-5);

        let pc_hadec = phase_centre.to_hadec(p.lmst);
        let ha_diff_arcmin = (pc_hadec.ha - p.hadec_j2000.ha).to_degrees() * 60.0;
        let dec_diff_arcmin = (pc_hadec.dec - p.hadec_j2000.dec).to_degrees() * 60.0;
        assert_abs_diff_eq!(ha_diff_arcmin, 11.034544, epsilon = 1e-4);
        assert_abs_diff_eq!(dec_diff_arcmin, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_precess_1099334672_to_j2000() {
        let epoch = Epoch::from_gpst_seconds(1099334672.0);
        let phase_centre = RADec::from_degrees(60.0, -30.0);

        let p = precess_time(
            MWA_LONG_RAD,
            MWA_LAT_RAD,
            phase_centre,
            epoch,
            Duration::from_total_nanoseconds(0),
        );
        let rm = precession_matrix(&p, phase_centre, epoch);
        // How do I know this is right? Good question! ... I don't.
        assert_abs_diff_eq!(rm[0][0], 1.0, epsilon = 1e-5);
        assert_abs_diff_eq!(rm[1][1], 1.0, epsilon = 1e-5);
        assert_abs_diff_eq!(rm[2][2], 1.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[0][1], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[0][2], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[1][0], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[1][2], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[2][0], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[2][1], 0.0, epsilon = 1e-5);

        // Expected J2000 values from ERFA: the zenith at the apparent sidereal
        // time (`eraGst06a`), taken to J2000 with `eraPnm06a`.
        assert_abs_diff_eq!(p.hadec_j2000.ha, 0.4100110770434626, epsilon = 1e-8);
        assert_abs_diff_eq!(p.hadec_j2000.dec, -0.5235987755982988, epsilon = 1e-10);
        assert_abs_diff_eq!(p.lmst, 1.4598017673520172, epsilon = 1e-10);
        assert_abs_diff_eq!(p.lmst_j2000, 1.4572086282400603, epsilon = 1e-8);
        assert_abs_diff_eq!(p.array_latitude_j2000, -0.46618075917013235, epsilon = 1e-8);
        assert_abs_diff_ne!(p.array_latitude_j2000, MWA_LAT_RAD, epsilon = 1e-5);

        let pc_hadec = phase_centre.to_hadec(p.lmst);
        let ha_diff_arcmin = (pc_hadec.ha - p.hadec_j2000.ha).to_degrees() * 60.0;
        let dec_diff_arcmin = (pc_hadec.dec - p.hadec_j2000.dec).to_degrees() * 60.0;
        assert_abs_diff_eq!(ha_diff_arcmin, 8.914556, epsilon = 1e-4);
        assert_abs_diff_eq!(dec_diff_arcmin, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_precess_1099334672_to_j2000_with_dut1() {
        // Test values have changed from the test above, due to the DUT1.

        let epoch = Epoch::from_gpst_seconds(1099334672.0);
        let phase_centre = RADec::from_degrees(60.0, -30.0);

        let p = precess_time(
            MWA_LONG_RAD,
            MWA_LAT_RAD,
            phase_centre,
            epoch,
            Duration::from_f64(-0.39623459, Unit::Second),
        );
        let rm = precession_matrix(&p, phase_centre, epoch);
        // How do I know this is right? Good question! ... I don't.
        assert_abs_diff_eq!(rm[0][0], 1.0, epsilon = 1e-5);
        assert_abs_diff_eq!(rm[1][1], 1.0, epsilon = 1e-5);
        assert_abs_diff_eq!(rm[2][2], 1.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[0][1], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[0][2], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[1][0], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[1][2], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[2][0], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(rm[2][1], 0.0, epsilon = 1e-5);

        // Expected J2000 values from ERFA: the zenith at the apparent sidereal
        // time (`eraGst06a`), taken to J2000 with `eraPnm06a`.
        assert_abs_diff_eq!(p.hadec_j2000.ha, 0.4099821814417699, epsilon = 1e-8);
        assert_abs_diff_eq!(p.hadec_j2000.dec, -0.5235987755982988, epsilon = 1e-10);
        assert_abs_diff_eq!(p.lmst, 1.4597728734531659, epsilon = 1e-10);
        assert_abs_diff_eq!(p.lmst_j2000, 1.4571797326383675, epsilon = 1e-8);
        assert_abs_diff_eq!(p.array_latitude_j2000, -0.4661808009369447, epsilon = 1e-8);
        assert_abs_diff_ne!(p.array_latitude_j2000, MWA_LAT_RAD, epsilon = 1e-5);

        let pc_hadec = phase_centre.to_hadec(p.lmst);
        let ha_diff_arcmin = (pc_hadec.ha - p.hadec_j2000.ha).to_degrees() * 60.0;
        let dec_diff_arcmin = (pc_hadec.dec - p.hadec_j2000.dec).to_degrees() * 60.0;
        assert_abs_diff_eq!(ha_diff_arcmin, 8.914561, epsilon = 1e-4);
        assert_abs_diff_eq!(dec_diff_arcmin, 0.0, epsilon = 1e-10);
    }
}
