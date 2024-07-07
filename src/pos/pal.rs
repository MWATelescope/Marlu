//! Useful Starlink PAL functions re-written in Rust.
//!
//! PAL depends on ERFA and (luckily) all the PAL functions we need just use
//! ERFA. Code is derived from <https://github.com/Starlink/pal> commit 7af65f0
//! and is licensed under the LGPL-3.0.
//!
//! Why not just depend on PAL? It's a huge pain. (1) The C library is either
//! named libstarlink-pal or libpal, (2) it depends on ERFA so statically
//! compiling it in a -sys crate is much harder than it should be, (3) this code
//! changes so slowly that we're unlikely to be out-of-date.

#![allow(non_snake_case)]
#![allow(clippy::excessive_precision)]

use erfa::{
    aliases::{
        eraAnp, eraC2s, eraEpj, eraEpj2jd, eraEpv00, eraGmst06, eraP06e, eraPdp, eraPmat06, eraPn,
        eraPnm06a, eraRx, eraRxp, eraRxpv, eraRxr, eraRz, eraS2c,
    },
    constants::{ERFA_AULT, ERFA_DAYSEC, ERFA_DJM0},
};

/// Greenwich mean sidereal time (consistent with IAU 2006 precession)
///
/// # Arguments
/// ut1 = double (Given)
///    Universal time (UT1) expressed as modified Julian Date (JD-2400000.5)
///
/// # Returned Value
/// Greenwich mean sidereal time
///
/// # Description
/// Greenwich mean sidereal time (consistent with IAU 2006 precession).
///
/// # Notes
/// - Uses `eraGmst06()`. See SOFA/ERFA documentation for details.
///
/// Original: <https://github.com/Starlink/pal/blob/7af65f05fcd33fd7362c586eae7e98972cb03f29/palOne2One.c#L1332>
pub fn palGmst(ut1: f64) -> f64 {
    eraGmst06(ERFA_DJM0, ut1, ERFA_DJM0, ut1)
}

/// Spherical coordinates to direction cosines
///
/// Arguments:
///    a = double (Given)
///       Spherical coordinate in radians (ra, long etc).
///    b = double (Given)
///       Spherical coordinate in radians (dec, lat etc).
///    v = double\[3\] (Returned)
///       x, y, z vector
/// Description:
///    The spherical coordinates are longitude (+ve anticlockwise looking
///    from the +ve latitude pole) and latitude.  The Cartesian coordinates
///    are right handed, with the x axis at zero longitude and latitude, and
///    the z axis at the +ve latitude pole.
/// Notes:
///    - Uses `eraS2c()`. See SOFA/ERFA documentation for details.
///
/// Original: <https://github.com/Starlink/pal/blob/7af65f05fcd33fd7362c586eae7e98972cb03f29/palOne2One.c#L368>
///
/// # Safety
///
/// v will be modified
pub unsafe fn palDcs2c(a: f64, b: f64, v: *mut f64) {
    let v2 = eraS2c(a, b);
    let v = std::slice::from_raw_parts_mut(v, 3);
    v.iter_mut().zip(v2).for_each(|(v, v2)| *v = v2);
}

/// Cartesian to spherical coordinates
///
/// Arguments:
///    v = double\[3\] (Given)
///       x, y, z vector.
///    a = double* (Returned)
///       Spherical coordinate (radians)
///    b = double* (Returned)
///       Spherical coordinate (radians)
/// Description:
///    The spherical coordinates are longitude (+ve anticlockwise looking
///    from the +ve latitude pole) and latitude.  The Cartesian coordinates
///    are right handed, with the x axis at zero longitude and latitude, and
///    the z axis at the +ve latitude pole.
/// Notes:
///    - Uses `eraC2s()`. See SOFA/ERFA documentation for details.
///
/// Original: <https://github.com/Starlink/pal/blob/7af65f05fcd33fd7362c586eae7e98972cb03f29/palOne2One.c#L327>
///
/// # Safety
///
/// `a` and `b` will be modified.
pub unsafe fn palDcc2s(v: *mut f64, a: &mut f64, b: &mut f64) {
    let v = std::slice::from_raw_parts(v, 3).try_into().unwrap();
    let (a2, b2) = eraC2s(v);
    *a = a2;
    *b = b2;
}

/// Normalize angle into range 0-2 pi
///
/// Arguments:
///    angle = double (Given)
///       angle in radians
/// Returned Value:
///    Angle expressed in the range 0-2 pi
/// Description:
///    Normalize angle into range 0-2 pi.
/// Notes:
///    - Uses `eraAnp()`. See SOFA/ERFA documentation for details.
///
/// Original: <https://github.com/Starlink/pal/blob/7af65f05fcd33fd7362c586eae7e98972cb03f29/palOne2One.c#L766>
pub fn palDranrm(angle: f64) -> f64 {
    eraAnp(angle)
}

/// Normalizes a 3-vector also giving the modulus
///
/// Arguments:
///    v = double\[3\] (Given)
///       vector
///    uv = double\[3\] (Returned)
///       unit vector in direction of "v"
///    vm = double* (Returned)
///       modulus of "v"
/// Description:
///    Normalizes a 3-vector also giving the modulus.
/// Notes:
///    - Uses `eraPn()`. See SOFA/ERFA documentation for details.
///    - the arguments are flipped
///
/// Original: <https://github.com/Starlink/pal/blob/7af65f05fcd33fd7362c586eae7e98972cb03f29/palOne2One.c#L1015>
///
/// # Safety
///
/// `uv` and `vm` will be modified.
pub unsafe fn palDvn(v: *mut f64, uv: *mut f64, vm: *mut f64) {
    let v = std::slice::from_raw_parts(v, 3).try_into().unwrap();
    let uv = std::slice::from_raw_parts_mut(uv, 3);
    let (m, uv2) = eraPn(v);
    *vm = m;
    uv.iter_mut().zip(uv2).for_each(|(uv, uv2)| *uv = uv2);
}

/// Scalar product of two 3-vectors
/// Arguments:
///    va = double\[3\] (Given)
///       First vector
///    vb = double\[3\] (Given)
///       Second vector
/// Returned Value:
///    Scalar product va.vb
/// Notes:
///    - Uses `eraPdp()`. See SOFA/ERFA documentation for details.
///
/// # Safety
///
/// [`erfa_sys::eraPdp`] requires `va` and `vb` to be mutable, even though
/// they are not mutated.
pub unsafe fn palDvdv(va: *mut f64, vb: *mut f64) -> f64 {
    let va: [f64; 3] = std::slice::from_raw_parts(va, 3).try_into().unwrap();
    let vb: [f64; 3] = std::slice::from_raw_parts(vb, 3).try_into().unwrap();
    eraPdp(va, vb)
}

/// Performs the 3-D forward unitary transformation
///
/// Arguments:
///    dm = double\[3\]\[3\] (Given)
///       matrix
///    va = double\[3\] (Given)
///       vector
///    dp = double\[3\] (Returned)
///       result vector
/// Notes:
///    - Uses `eraRxp()`. See SOFA/ERFA documentation for details.
///
/// # Safety
///
/// `dp` will be modified.
pub unsafe fn palDmxv(dm: *mut [f64; 3], va: *mut f64, vb: *mut f64) {
    let dm: [[f64; 3]; 3] = std::slice::from_raw_parts(dm, 3).try_into().unwrap();
    let va: [f64; 3] = std::slice::from_raw_parts(va, 3).try_into().unwrap();
    let vb = std::slice::from_raw_parts_mut(vb, 3);
    let vb2 = eraRxp(dm, va);
    vb.iter_mut().zip(vb2).for_each(|(vb, vb2)| *vb = vb2);
}

/// Barycentric and heliocentric velocity and position of the Earth.
///
/// Arguments:
///    date = double (Given)
///       TDB (loosely ET) as a Modified Julian Date (JD-2400000.5)
///    deqx = double (Given)
///       Julian epoch (e.g. 2000.0) of mean equator and equinox of the
///       vectors returned.  If deqx <= 0.0, all vectors are referred to the
///       mean equator and equinox (FK5) of epoch date.
///    dvb = double\[3\] (Returned)
///       Barycentric velocity (AU/s, AU)
///    dpb = double\[3\] (Returned)
///       Barycentric position (AU/s, AU)
///    dvh = double\[3\] (Returned)
///       heliocentric velocity (AU/s, AU)
///    dph = double\[3\] (Returned)
///       Heliocentric position (AU/s, AU)
/// Description:
///    Returns the barycentric and heliocentric velocity and position of the
///    Earth at a given epoch, given with respect to a specified equinox.
///    For information about accuracy, see the function eraEpv00.
///
/// # Safety
///
/// `dvb`, `dpb`, `dvh`, `dph` will be modified.
pub unsafe fn palEvp(
    date: f64,
    deqx: f64,
    dvb: *mut f64,
    dpb: *mut f64,
    dvh: *mut f64,
    dph: *mut f64,
) {
    /* BCRS PV-vectors. */
    // CHJ: PAL suppresses a warning indicated by eraEpv00, if the date is not
    // within 1900-2100.
    let (_, mut pvh, mut pvb) = eraEpv00(2400000.5, date);

    /* Was precession to another equinox requested? */
    if deqx > 0.0 {
        /* Yes: compute precession matrix from J2000.0 to deqx. */
        let (d1, d2) = eraEpj2jd(deqx);
        let r = eraPmat06(d1, d2);

        /* Rotate the PV-vectors. */
        let pvh2 = eraRxpv(r, pvh);
        pvh = pvh2;
        let pvb2 = eraRxpv(r, pvb);
        pvb = pvb2;
    }

    /* Return the required vectors. */
    let dvb = std::slice::from_raw_parts_mut(dvb, 3);
    let dpb = std::slice::from_raw_parts_mut(dpb, 3);
    let dvh = std::slice::from_raw_parts_mut(dvh, 3);
    let dph = std::slice::from_raw_parts_mut(dph, 3);
    for i in 0..3 {
        dvh[i] = pvh[1][i] / ERFA_DAYSEC;
        dvb[i] = pvb[1][i] / ERFA_DAYSEC;
        dph[i] = pvh[0][i];
        dpb[i] = pvb[0][i];
    }
}

/// Form the matrix of bias-precession-nutation (IAU 2006/2000A)
///
///
/// # Arguments
///
/// epoch = double (Returned)
///    Julian epoch for mean coordinates.
/// date = double (Returned)
///    Modified Julian Date (JD-2400000.5) for true coordinates.
/// rmatpn = double\[3\]\[3\] (Returned)
///    combined NPB matrix
///
/// # Description
///
/// Form the matrix of bias-precession-nutation (IAU 2006/2000A).
/// The epoch and date are TT (but TDB is usually close enough).
/// The matrix is in the sense   v(true)  =  rmatpn * v(mean).
///
/// Original: <https://github.com/Starlink/pal/blob/master/palPrenut.c>
///
/// # Safety
///
/// `epoch`, `date`, and `rmatpn` will be modified.
pub unsafe fn palPrenut(epoch: f64, date: f64, rmatpn: *mut [f64; 3]) {
    /* Specified Julian epoch as a 2-part JD. */
    let (d1, d2) = eraEpj2jd(epoch);

    /* P matrix, from specified epoch to J2000.0. */
    let (eps0, psia, oma, _, _, _, _, _, chia, _, _, _, _, _, _, _) = eraP06e(d1, d2);
    let mut r1 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    eraRz(-chia, &mut r1);
    eraRx(oma, &mut r1);
    eraRz(psia, &mut r1);
    eraRx(-eps0, &mut r1);

    /* NPB matrix, from J2000.0 to date. */
    let r2 = eraPnm06a(ERFA_DJM0, date);

    /* NPB matrix, from specified epoch to date. */
    let rmatpn2 = eraRxr(r2, r1);
    let rmatpn = std::slice::from_raw_parts_mut(rmatpn, 3);
    rmatpn[0][0] = rmatpn2[0][0];
    rmatpn[0][1] = rmatpn2[0][1];
    rmatpn[0][2] = rmatpn2[0][2];
    rmatpn[1][0] = rmatpn2[1][0];
    rmatpn[1][1] = rmatpn2[1][1];
    rmatpn[1][2] = rmatpn2[1][2];
    rmatpn[2][0] = rmatpn2[2][0];
    rmatpn[2][1] = rmatpn2[2][1];
    rmatpn[2][2] = rmatpn2[2][2];
}

/// Compute parameters needed by palAmpqk and palMapqk.
///
/// Arguments:
///    eq = double (Given)
///       epoch of mean equinox to be used (Julian)
///    date = double (Given)
///       TDB (JD-2400000.5)
///    amprms =   double\[21\]  (Returned)
///       star-independent mean-to-apparent parameters:
///       - (0)      time interval for proper motion (Julian years)
///       - (1-3)    barycentric position of the Earth (AU)
///       - (4-6)    heliocentric direction of the Earth (unit vector)
///       - (7)      (Schwarzschild radius of Sun)/(Sun-Earth distance)
///       - (8-10)   abv: barycentric Earth velocity in units of c
///       - (11)     sqrt(1-v^2) where v=modulus(abv)
///       - (12-20)  precession/nutation (3,3) matrix
///
/// Description:
///    Compute star-independent parameters in preparation for
///    transformations between mean place and geocentric apparent place.
///    The parameters produced by this function are required in the
///    parallax, aberration, and nutation/bias/precession parts of the
///    mean/apparent transformations.
///    The reference systems and timescales used are IAU 2006.
///
/// # Safety
///
/// `amprms` will be modified.
pub unsafe fn palMappa(eq: f64, date: f64, amprms: *mut f64) {
    /* Local constants */

    /*  Gravitational radius of the Sun x 2 (2*mu/c**2, AU) */
    let gr2: f64 = 2.0 * 9.87063e-9;

    /* Local Variables; */
    let mut ebd = [0.0; 3];
    let mut ehd = [0.0; 3];
    let mut eh = [0.0; 3];

    /* Initialise so that unsused values are returned holding zero */
    let amprms = std::slice::from_raw_parts_mut(amprms, 21);
    amprms.fill(0.0);

    /* Time interval for proper motion correction. */
    amprms[0] = eraEpj(ERFA_DJM0, date) - eq;

    /* Get Earth barycentric and heliocentric position and velocity. */
    palEvp(
        date,
        eq,
        ebd.as_mut_ptr(),
        &mut amprms[1],
        ehd.as_mut_ptr(),
        eh.as_mut_ptr(),
    );

    /* Heliocentric direction of Earth (normalized) and modulus. */
    let (e, amprms_from_4) = eraPn(eh);
    amprms[4] = amprms_from_4[0];
    amprms[5] = amprms_from_4[1];
    amprms[6] = amprms_from_4[2];

    /* Light deflection parameter */
    amprms[7] = gr2 / e;

    /* Aberration parameters. */
    for i in 0..3 {
        amprms[i + 8] = ebd[i] * ERFA_AULT;
    }
    let (vm, _) = eraPn(amprms[8..11].try_into().unwrap());
    amprms[11] = (1.0 - vm * vm).sqrt();

    /* NPB matrix. */
    palPrenut(eq, date, amprms.as_mut_ptr().add(12) as _);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ndarray::prelude::*;
    use approx::assert_abs_diff_eq;

    /// We use the SOFA/ERFA test values rather than the values from SLA
    /// because the precession models have changed
    ///
    /// Original: <https://github.com/Starlink/pal/blob/7af65f05fcd33fd7362c586eae7e98972cb03f29/palTest.c#L941>
    #[test]
    fn gmst() {
        assert_abs_diff_eq!(palGmst(53736.), 1.754174971870091203, epsilon = 1e-12);
    }

    /// Test all the 3-vector and 3x3 matrix routines.
    ///
    /// Original: <https://github.com/Starlink/pal/blob/7af65f05fcd33fd7362c586eae7e98972cb03f29/palTest.c#L1148>
    #[test]
    fn vecmat() {
        // let mut dav = [-0.123, 0.0987, 0.0654];

        //     /* palDcs2c et al */
        let mut dv1 = [0_f64; 3];
        let mut dv2 = [0_f64; 3];
        let mut dv3 = [0_f64; 3];
        // let mut dv4 = [ -0.5366267667260526, 0.06977111097651445, -0.8409302618566215 ];
        let mut dv5 = [0.006889040510209034, -1.577473205461961, 0.5201843672856759];
        let mut dv6 = [0_f64; 3];
        // let mut dv7 = [0_f64; 3];
        let mut dvm = 0_f64;

        // let mut drm: [[f64; 3]; 3] = [
        //     [-0.09010460088585805, 0.3075993402463796, 0.9472400998581048],
        //     [-0.3161868071070688, 0.8930686362478707, -0.3200848543149236],
        //     [ -0.9444083141897035, -0.3283459407855694, 0.01678926022795169],
        // ];

        let mut drm1: [[f64; 3]; 3] = [
            [0.9930075842721269, 0.05902743090199868, -0.1022335560329612],
            [
                -0.07113807138648245,
                0.9903204657727545,
                -0.1191836812279541,
            ],
            [0.09420887631983825, 0.1256229973879967, 0.9875948309655174],
        ];

        let mut drm2: [[f64; 3]; 3] = [
            [-0.1681574770810878, 0.1981362273264315, 0.9656423242187410],
            [-0.2285369373983370, 0.9450659587140423, -0.2337117924378156],
            [
                -0.9589024617479674,
                -0.2599853247796050,
                -0.1136384607117296,
            ],
        ];

        /* palDcs2c */
        unsafe { palDcs2c(3.0123, -0.999, dv1.as_mut_ptr()) };
        assert_abs_diff_eq!(dv1[0], -0.5366267667260525, epsilon = 1e-12);
        assert_abs_diff_eq!(dv1[1], 0.06977111097651444, epsilon = 1e-12);
        assert_abs_diff_eq!(dv1[2], -0.8409302618566215, epsilon = 1e-12);

        /* palDmxv */
        unsafe { palDmxv(drm1.as_mut_ptr(), dv1.as_mut_ptr(), dv2.as_mut_ptr()) };
        unsafe { palDmxv(drm2.as_mut_ptr(), dv2.as_mut_ptr(), dv3.as_mut_ptr()) };
        assert_abs_diff_eq!(dv3[0], -0.7267487768696160, epsilon = 1e-10);
        assert_abs_diff_eq!(dv3[1], 0.5011537352639822, epsilon = 1e-12);
        assert_abs_diff_eq!(dv3[2], 0.4697671220397141, epsilon = 1e-12);

        dv5.iter_mut().for_each(|x| *x *= 1000.0);

        /* palDvn */
        unsafe { palDvn(dv5.as_mut_ptr(), dv6.as_mut_ptr(), &mut dvm) };
        assert_abs_diff_eq!(dv6[0], 0.004147420704640065, epsilon = 1e-12);
        assert_abs_diff_eq!(dv6[1], -0.9496888606842218, epsilon = 1e-12);
        assert_abs_diff_eq!(dv6[2], 0.3131674740355448, epsilon = 1e-12);
        assert_abs_diff_eq!(dvm, 1661.042127339937, epsilon = 1e-9);
    }

    /// Test palDcc2s
    ///
    /// Original <https://github.com/Starlink/pal/blob/7af65f05fcd33fd7362c586eae7e98972cb03f29/palTest.c#L474>
    #[test]
    fn dcc2s() {
        let mut dv = [100_f64, -50_f64, 25_f64];
        let mut da = 0_f64;
        let mut db = 0_f64;

        unsafe { palDcc2s(dv.as_mut_ptr(), &mut da, &mut db) };
        assert_abs_diff_eq!(da, -0.4636476090008061, epsilon = 1e-12);
        assert_abs_diff_eq!(db, 0.2199879773954594, epsilon = 1e-12);
    }

    /// Original: <https://github.com/Starlink/pal/blob/7af65f05fcd33fd7362c586eae7e98972cb03f29/palTest.c#L1074>
    #[test]
    fn dranrm() {
        assert_abs_diff_eq!(palDranrm(-0.1), 6.183185307179587, epsilon = 1e-12);
    }

    /// Original: <https://github.com/Starlink/pal/blob/7af65f05fcd33fd7362c586eae7e98972cb03f29/palTest.c#L1329>
    #[test]
    fn evp() {
        let mut dvb = Array1::zeros(3);
        let mut dpb = Array1::zeros(3);
        let mut dvh = Array1::zeros(3);
        let mut dph = Array1::zeros(3);

        let vbex = array![
            1.6957348127008098514e-07,
            -9.1093446116039685966e-08,
            -3.9528532243991863036e-08,
        ];
        let pbex = array![
            -0.49771075259730546136,
            -0.80273812396332311359,
            -0.34851593942866060383,
        ];
        let vhex = array![
            1.6964379181455713805e-07,
            -9.1147224045727438391e-08,
            -3.9553158272334222497e-08,
        ];
        let phex = array![
            -0.50169124421419830639,
            -0.80650980174901798492,
            -0.34997162028527262212,
        ];

        unsafe {
            palEvp(
                2010.0,
                2012.0,
                dvb.as_mut_ptr(),
                dpb.as_mut_ptr(),
                dvh.as_mut_ptr(),
                dph.as_mut_ptr(),
            );
        };

        assert_abs_diff_eq!(dvb, vbex, epsilon = 1e-12);
        assert_abs_diff_eq!(dpb, pbex, epsilon = 1e-12);
        assert_abs_diff_eq!(dvh, vhex, epsilon = 1e-12);
        assert_abs_diff_eq!(dph, phex, epsilon = 1e-12);
    }

    #[test]
    fn prenut() {
        // TODO: PAL doesn't test this.
    }

    /// Original: <https://github.com/Starlink/pal/blob/7af65f05fcd33fd7362c586eae7e98972cb03f29/palTest.c#L1395>
    #[test]
    fn mappa() {
        let mut amprms = Array1::zeros(21);
        let expected = array![
            1.9986310746064646082_f64,
            -0.1728200754134739392,
            0.88745394651412767839,
            0.38472374350184274094,
            -0.17245634725219796679,
            0.90374808622520386159,
            0.3917884696321610738,
            2.0075929387510784968e-08,
            -9.9464149073251757597e-05,
            -1.6125306981057062306e-05,
            -6.9897255793245634435e-06,
            0.99999999489900059935,
            0.99999983777998024959,
            -0.00052248206600935195865,
            -0.00022683144398381763045,
            0.00052248547063364874764,
            0.99999986339269864022,
            1.4950491424992534218e-05,
            0.00022682360163333854623,
            -1.5069005133483779417e-05,
            0.99999997416198904698,
        ];

        unsafe { palMappa(2010.0, 55927.0, amprms.as_mut_ptr()) };

        assert_abs_diff_eq!(amprms, expected, epsilon = 1e-12);
    }
}
