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

use hifitime::{Duration, Epoch};

use crate::{pal, HADec, RADec, XyzGeodetic};

#[derive(Debug)]
pub struct PrecessionInfo {
    /// Bias procession rotation matrix.
    rotation_matrix: [[f64; 3]; 3],

    /// The precessed phase centre in the J2000 epoch.
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
/// This function is very similar to cotter's `PrepareTimestepUVW`.
pub fn precess_time(
    array_longitude_rad: f64,
    array_latitude_rad: f64,
    phase_centre: RADec,
    time: Epoch,
    dut1: Duration,
) -> PrecessionInfo {
    // Note that we explicitly use the mean LST (i.e. LMST) because we're
    // handling nutation ourselves.
    let lmst = get_lmst(array_longitude_rad, time, dut1);

    let j2000 = 2000.0;
    let mjd = (time + dut1).to_mjd_utc_days();
    let radec_aber = aber_radec_rad(j2000, mjd, phase_centre);
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

    let precessed = hadec_j2000(&mut rotation_matrix, lmst, array_latitude_rad, radec_aber);

    PrecessionInfo {
        rotation_matrix,
        hadec_j2000: precessed.hadec,
        lmst,
        lmst_j2000: precessed.lmst,
        array_latitude_j2000: precessed.latitude,
    }
}

// Blatently stolen from cotter.
fn aber_radec_rad(eq: f64, mjd: f64, radec: RADec) -> RADec {
    let mut v1 = [0.0; 3];
    let mut v2 = [0.0; 3];

    unsafe {
        pal::palDcs2c(radec.ra, radec.dec, v1.as_mut_ptr());
        stelaber(eq, mjd, &mut v1, &mut v2);
        let mut ra2 = 0.0;
        let mut dec2 = 0.0;
        pal::palDcc2s(v2.as_mut_ptr(), &mut ra2, &mut dec2);
        ra2 = pal::palDranrm(ra2);

        RADec::from_radians(ra2, dec2)
    }
}

// Blatently stolen from cotter.
fn stelaber(eq: f64, mjd: f64, v1: &mut [f64; 3], v2: &mut [f64; 3]) {
    let mut amprms = [0.0; 21];
    let mut v1n = [0.0; 3];
    let mut v2un = [0.0; 3];
    let mut abv = [0.0; 3];

    unsafe {
        pal::palMappa(eq, mjd, amprms.as_mut_ptr());

        /* Unpack scalar and vector parameters */
        let ab1 = &amprms[11];
        abv[0] = amprms[8];
        abv[1] = amprms[9];
        abv[2] = amprms[10];

        let mut w = 0.0;
        pal::palDvn(v1.as_mut_ptr(), v1n.as_mut_ptr(), &mut w);

        /* Aberration (normalization omitted) */
        let p1dv = pal::palDvdv(v1n.as_mut_ptr(), abv.as_mut_ptr());
        w = 1.0 + p1dv / (ab1 + 1.0);
        v2un[0] = ab1 * v1n[0] + w * abv[0];
        v2un[1] = ab1 * v1n[1] + w * abv[1];
        v2un[2] = ab1 * v1n[2] + w * abv[2];

        /* Normalize */
        pal::palDvn(v2un.as_mut_ptr(), v2.as_mut_ptr(), &mut w);
    }
}

/// The return type of `hadec_j2000`. All values are in radians.
struct HADecJ2000 {
    hadec: HADec,
    latitude: f64,
    lmst: f64,
}

fn hadec_j2000(
    rotation_matrix: &mut [[f64; 3]; 3],
    lmst: f64,
    lat_rad: f64,
    radec: RADec,
) -> HADecJ2000 {
    let (new_lmst, new_lat) = rotate_radec(rotation_matrix, lmst, lat_rad);
    HADecJ2000 {
        hadec: HADec::from_radians(pal::palDranrm(new_lmst - radec.ra), radec.dec),
        latitude: new_lat,
        lmst: new_lmst,
    }
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
        assert_abs_diff_eq!(eor0.rotation_matrix[0][0], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(eor0.rotation_matrix[1][1], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(eor0.rotation_matrix[2][2], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(eor0.rotation_matrix[0][1], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(eor0.rotation_matrix[0][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(eor0.rotation_matrix[1][0], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(eor0.rotation_matrix[1][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(eor0.rotation_matrix[2][0], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(eor0.rotation_matrix[2][1], 0.0, epsilon = 3e-5);

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
        assert_abs_diff_eq!(eor1.rotation_matrix[0][0], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(eor1.rotation_matrix[1][1], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(eor1.rotation_matrix[2][2], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(eor1.rotation_matrix[0][1], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(eor1.rotation_matrix[0][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(eor1.rotation_matrix[1][0], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(eor1.rotation_matrix[1][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(eor1.rotation_matrix[2][0], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(eor1.rotation_matrix[2][1], 0.0, epsilon = 3e-5);

        assert_abs_diff_eq!(
            eor1.lmst_j2000 - eor1.hadec_j2000.ha,
            -5.235898085317921,
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
        assert_abs_diff_eq!(eor0.rotation_matrix[0][0], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(eor0.rotation_matrix[1][1], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(eor0.rotation_matrix[2][2], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(eor0.rotation_matrix[0][1], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(eor0.rotation_matrix[0][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(eor0.rotation_matrix[1][0], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(eor0.rotation_matrix[1][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(eor0.rotation_matrix[2][0], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(eor0.rotation_matrix[2][1], 0.0, epsilon = 3e-5);

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
        assert_abs_diff_eq!(eor1.rotation_matrix[0][0], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(eor1.rotation_matrix[1][1], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(eor1.rotation_matrix[2][2], 1.0, epsilon = 1e-8);
        assert_abs_diff_eq!(eor1.rotation_matrix[0][1], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(eor1.rotation_matrix[0][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(eor1.rotation_matrix[1][0], 0.0, epsilon = 7e-5);
        assert_abs_diff_eq!(eor1.rotation_matrix[1][2], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(eor1.rotation_matrix[2][0], 0.0, epsilon = 3e-5);
        assert_abs_diff_eq!(eor1.rotation_matrix[2][1], 0.0, epsilon = 3e-5);

        assert_abs_diff_eq!(
            eor1.lmst_j2000 - eor1.hadec_j2000.ha,
            -5.235898085317921,
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
        // How do I know this is right? Good question! ... I don't.
        assert_abs_diff_eq!(p.rotation_matrix[0][0], 1.0, epsilon = 1e-5);
        assert_abs_diff_eq!(p.rotation_matrix[1][1], 1.0, epsilon = 1e-5);
        assert_abs_diff_eq!(p.rotation_matrix[2][2], 1.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[0][1], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[0][2], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[1][0], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[1][2], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[2][0], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[2][1], 0.0, epsilon = 1e-5);

        assert_abs_diff_eq!(p.hadec_j2000.ha, 6.0714305189419715, epsilon = 1e-10);
        assert_abs_diff_eq!(p.hadec_j2000.dec, -0.47122418312765446, epsilon = 1e-10);
        assert_abs_diff_eq!(p.lmst, 6.0747789094260245, epsilon = 1e-10);
        assert_abs_diff_eq!(p.lmst_j2000, 6.071524853456497, epsilon = 1e-10);
        assert_abs_diff_eq!(p.array_latitude_j2000, -0.467396549790915, epsilon = 1e-10);
        assert_abs_diff_ne!(p.array_latitude_j2000, MWA_LAT_RAD, epsilon = 1e-5);

        let pc_hadec = phase_centre.to_hadec(p.lmst);
        let ha_diff_arcmin = (pc_hadec.ha - p.hadec_j2000.ha).to_degrees() * 60.0;
        let dec_diff_arcmin = (pc_hadec.dec - p.hadec_j2000.dec).to_degrees() * 60.0;
        assert_abs_diff_eq!(ha_diff_arcmin, 11.510918573880216, epsilon = 1e-5); // 11.5 arcmin!
        assert_abs_diff_eq!(dec_diff_arcmin, -0.05058613713495692, epsilon = 1e-5);
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
        // How do I know this is right? Good question! ... I don't.
        assert_abs_diff_eq!(p.rotation_matrix[0][0], 1.0, epsilon = 1e-5);
        assert_abs_diff_eq!(p.rotation_matrix[1][1], 1.0, epsilon = 1e-5);
        assert_abs_diff_eq!(p.rotation_matrix[2][2], 1.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[0][1], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[0][2], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[1][0], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[1][2], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[2][0], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[2][1], 0.0, epsilon = 1e-5);

        assert_abs_diff_eq!(p.hadec_j2000.ha, 0.409885996082088, epsilon = 1e-10);
        assert_abs_diff_eq!(p.hadec_j2000.dec, -0.5235637661235192, epsilon = 1e-10);
        assert_abs_diff_eq!(p.lmst, 1.4598017673520172, epsilon = 1e-10);
        assert_abs_diff_eq!(p.lmst_j2000, 1.4571918352968762, epsilon = 1e-10);
        assert_abs_diff_eq!(p.array_latitude_j2000, -0.4661807836570052, epsilon = 1e-10);
        assert_abs_diff_ne!(p.array_latitude_j2000, MWA_LAT_RAD, epsilon = 1e-5);

        let pc_hadec = phase_centre.to_hadec(p.lmst);
        let ha_diff_arcmin = (pc_hadec.ha - p.hadec_j2000.ha).to_degrees() * 60.0;
        let dec_diff_arcmin = (pc_hadec.dec - p.hadec_j2000.dec).to_degrees() * 60.0;
        assert_abs_diff_eq!(ha_diff_arcmin, 9.344552279378359, epsilon = 1e-5);
        assert_abs_diff_eq!(dec_diff_arcmin, -0.12035370887056628, epsilon = 1e-5);
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
        // How do I know this is right? Good question! ... I don't.
        assert_abs_diff_eq!(p.rotation_matrix[0][0], 1.0, epsilon = 1e-5);
        assert_abs_diff_eq!(p.rotation_matrix[1][1], 1.0, epsilon = 1e-5);
        assert_abs_diff_eq!(p.rotation_matrix[2][2], 1.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[0][1], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[0][2], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[1][0], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[1][2], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[2][0], 0.0, epsilon = 1e-5);
        assert_abs_diff_ne!(p.rotation_matrix[2][1], 0.0, epsilon = 1e-5);

        assert_abs_diff_eq!(p.hadec_j2000.ha, 0.4098571004815674, epsilon = 1e-10);
        assert_abs_diff_eq!(p.hadec_j2000.dec, -0.5235637661235192, epsilon = 1e-10);
        assert_abs_diff_eq!(p.lmst, 1.4597728734531659, epsilon = 1e-10);
        assert_abs_diff_eq!(p.lmst_j2000, 1.457162939693978, epsilon = 1e-10);
        assert_abs_diff_eq!(p.array_latitude_j2000, -0.4661808254250147, epsilon = 1e-10);
        assert_abs_diff_ne!(p.array_latitude_j2000, MWA_LAT_RAD, epsilon = 1e-5);

        let pc_hadec = phase_centre.to_hadec(p.lmst);
        let ha_diff_arcmin = (pc_hadec.ha - p.hadec_j2000.ha).to_degrees() * 60.0;
        let dec_diff_arcmin = (pc_hadec.dec - p.hadec_j2000.dec).to_degrees() * 60.0;
        assert_abs_diff_eq!(ha_diff_arcmin, 9.344552279378359, epsilon = 1e-5);
        assert_abs_diff_eq!(dec_diff_arcmin, -0.12035370887056628, epsilon = 1e-5);
    }
}
