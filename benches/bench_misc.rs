// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Misc Benchmarks

use criterion::*;
use marlu::{
    c64,
    ndarray::{Array1, Array3},
    pos::xyz,
    HADec, Jones, XyzGeodetic,
};

// /////////////////////// //
// Miscelaneous Benchmarks //
// /////////////////////// //

#[inline]
fn mul(j1: [c64; 4], j2: [c64; 4]) -> [c64; 4] {
    [
        j1[0] * j2[0] + j1[1] * j2[2],
        j1[0] * j2[1] + j1[1] * j2[3],
        j1[2] * j2[0] + j1[3] * j2[2],
        j1[2] * j2[1] + j1[3] * j2[3],
    ]
}

#[inline]
fn mul2(j1: &[c64; 4], j2: &[c64; 4]) -> [c64; 4] {
    [
        j1[0] * j2[0] + j1[1] * j2[2],
        j1[0] * j2[1] + j1[1] * j2[3],
        j1[2] * j2[0] + j1[3] * j2[2],
        j1[2] * j2[1] + j1[3] * j2[3],
    ]
}

fn misc(c: &mut Criterion) {
    // Is the parallel xyzs_to_uvws really worth it?
    c.bench_function("xyzs_to_uvws", |b| {
        // The values are irrelevant.
        let xyzs = vec![XyzGeodetic::default(); 8128];
        let phase_centre = HADec::new(0.0, -27.0);
        b.iter(|| xyz::xyzs_to_uvws(&xyzs, phase_centre))
    });

    c.bench_function("xyzs_to_uvws_parallel", |b| {
        // The values are irrelevant.
        let xyzs = vec![XyzGeodetic::default(); 8128];
        let phase_centre = HADec::new(0.0, -27.0);
        b.iter(|| xyz::xyzs_to_uvws_parallel(&xyzs, phase_centre))
    });

    // Sanity check that allocating Jones<f64> has no overhead compared to [f64;
    // 8].
    c.bench_function("allocating many Jones<f64>", |b| {
        b.iter(|| {
            let shape = (20, 8128, 384);
            let _a: Array3<Jones<f64>> = Array3::from_elem(shape, Jones::default());
        })
    });

    c.bench_function("allocating many [c64; 4]", |b| {
        b.iter(|| {
            let shape = (20, 8128, 384);
            let _a: Array3<[c64; 4]> = Array3::from_elem(shape, [c64::new(0.0, 0.0); 4]);
        })
    });

    // Sanity check that multiplying Jones<f64> has no overhead compared to
    // multiplying [c64; 4].
    c.bench_function("multiply Jones<f64>", |b| {
        let i = c64::new(1.0, 2.0);
        let j1 = Jones::from([i, i + 1.0, i + 2.0, i + 3.0]);
        let j2 = Jones::from([i * 2.0, i * 3.0, i * 4.0, i * 5.0]);
        b.iter(|| {
            let _j3 = black_box(j1 * j2);
        })
    });

    c.bench_function("multiply [c64; 4]", |b| {
        let i = c64::new(1.0, 2.0);
        let j1 = [i, i + 1.0, i + 2.0, i + 3.0];
        let j2 = [i * 2.0, i * 3.0, i * 4.0, i * 5.0];
        b.iter(|| {
            let _j3 = black_box(mul(j1, j2));
        })
    });

    c.bench_function("multiply Array1<Jones<f64>>", |b| {
        let i = c64::new(1.0, 2.0);
        let a1 = Array1::from_elem(1000000, Jones::from([i, i + 1.0, i + 2.0, i + 3.0]));
        let a2 = Array1::from_elem(1000000, Jones::from([i * 2.0, i * 3.0, i * 4.0, i * 5.0]));
        b.iter(|| {
            let _a3 = black_box(&a1 * &a2);
        })
    });

    c.bench_function("multiply Array1<[c64; 4]>", |b| {
        let i = c64::new(1.0, 2.0);
        let a1 = Array1::from_elem(1000000, [i, i + 1.0, i + 2.0, i + 3.0]);
        let a2 = Array1::from_elem(1000000, [i * 2.0, i * 3.0, i * 4.0, i * 5.0]);
        b.iter(|| {
            let _a3: Array1<[c64; 4]> = a1
                .iter()
                .zip(a2.iter())
                .map(|(a1, a2)| black_box(mul2(a1, a2)))
                .collect();
        })
    });
}

criterion_group!(benches, misc);
criterion_main!(benches);
