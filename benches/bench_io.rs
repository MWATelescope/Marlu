// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Benchmarks

use criterion::*;
use glob::glob;
use marlu::{
    io::{ms::MeasurementSetWriter, VisWritable},
    mwalib::CorrelatorContext,
    ndarray::Array3,
    Complex, Jones, RADec,
};
use mwalib::get_baseline_count;
use ndarray::Array4;
use std::env;
use std::path::Path;
use tempfile::tempdir;

// ///////////// //
// IO Benchmarks //
// ///////////// //

const NUM_TIMESTEPS: usize = 2;
const NUM_COARSE_CHANS: usize = 2;
const NUM_ITERATIONS: usize = 10;

fn get_test_dir() -> String {
    env::var("MARLU_TEST_DIR").unwrap_or_else(|_| String::from("/mnt/data"))
}

fn get_context_mwax_half_1247842824() -> CorrelatorContext {
    let test_dir = get_test_dir();
    let test_path = Path::new(&test_dir);
    let vis_path = test_path.join("1247842824_vis");
    let metafits_path = vis_path
        .join("1247842824.metafits")
        .to_str()
        .unwrap()
        .to_owned();
    let gpufits_glob = vis_path
        .join("1247842824_*gpubox*_00.fits")
        .to_str()
        .unwrap()
        .to_owned();
    let gpufits_files: Vec<String> = glob(gpufits_glob.as_str())
        .unwrap()
        .filter_map(Result::ok)
        .map(|path| path.to_str().unwrap().to_owned())
        .collect();
    CorrelatorContext::new(&metafits_path, &gpufits_files).unwrap()
}

fn bench_ms_init_mwax_half_1247842824(crt: &mut Criterion) {
    let context = get_context_mwax_half_1247842824();
    let first_common_timestep_idx = *(context.common_timestep_indices.first().unwrap());
    let mwalib_timestep_range =
        first_common_timestep_idx..(first_common_timestep_idx + NUM_TIMESTEPS);
    let mwalib_coarse_chan_range = 0..NUM_COARSE_CHANS;

    let tmp_dir = tempdir().unwrap();
    let ms_path = tmp_dir.path().join("1254670392.none.ms");

    let phase_centre = RADec::from_mwalib_phase_or_pointing(&context.metafits_context);
    let ms_writer = MeasurementSetWriter::new(ms_path, phase_centre, None);

    crt.bench_function(
        "ms_writer::initialize_from_mwalib - mwax_half_1247842824",
        |bch| {
            bch.iter(|| {
                ms_writer
                    .initialize_from_mwalib(
                        &context,
                        &mwalib_timestep_range,
                        &mwalib_coarse_chan_range,
                    )
                    .unwrap();
            })
        },
    );
}

fn synthesize_test_data(
    shape: (usize, usize, usize),
) -> (Array3<Jones<f32>>, Array4<f32>, Array4<bool>) {
    let jones_array = Array3::from_shape_fn(shape, |(timestep_idx, chan_idx, baseline_idx)| {
        Jones::from([
            Complex::new(0., timestep_idx as _),
            Complex::new(0., chan_idx as _),
            Complex::new(0., baseline_idx as _),
            Complex::new(0., 1.),
        ])
    });
    let shape_with_pol = (shape.0, shape.1, shape.2, 4);
    let weight_array = Array4::from_shape_fn(
        shape_with_pol,
        |(timestep_idx, chan_idx, baseline_idx, pol_idx)| {
            (timestep_idx
                * shape_with_pol.0
                * shape_with_pol.1
                * shape_with_pol.2
                * shape_with_pol.3
                + chan_idx * shape_with_pol.0 * shape_with_pol.1 * shape_with_pol.2
                + baseline_idx * shape_with_pol.0 * shape_with_pol.1
                + pol_idx) as f32
        },
    );
    let flag_array = Array4::from_shape_fn(
        shape_with_pol,
        |(timestep_idx, chan_idx, baseline_idx, pol_idx)| {
            (timestep_idx
                * shape_with_pol.0
                * shape_with_pol.1
                * shape_with_pol.2
                * shape_with_pol.3
                + chan_idx * shape_with_pol.0 * shape_with_pol.1 * shape_with_pol.2
                + baseline_idx * shape_with_pol.0 * shape_with_pol.1
                + pol_idx % 2
                == 0) as _
        },
    );
    (jones_array, weight_array, flag_array)
}

fn bench_ms_write_mwax_part_1247842824(crt: &mut Criterion) {
    let context = get_context_mwax_half_1247842824();
    let first_common_timestep_idx = *(context.common_timestep_indices.first().unwrap());
    let mwalib_timestep_range =
        first_common_timestep_idx..(first_common_timestep_idx + NUM_TIMESTEPS);
    let mwalib_coarse_chan_range = 0..NUM_COARSE_CHANS;

    let tmp_dir = tempdir().unwrap();
    let ms_path = tmp_dir.path().join("1254670392.none.ms");

    let phase_centre = RADec::from_mwalib_phase_or_pointing(&context.metafits_context);
    let mut ms_writer = MeasurementSetWriter::new(ms_path, phase_centre, None);

    ms_writer
        .initialize_from_mwalib(&context, &mwalib_timestep_range, &mwalib_coarse_chan_range)
        .unwrap();

    let mwalib_baseline_idxs: Vec<usize> =
        (0..get_baseline_count(context.metafits_context.num_ants)).collect();

    let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;

    let shape = (
        mwalib_timestep_range.len(),
        mwalib_coarse_chan_range.len() * fine_chans_per_coarse,
        mwalib_baseline_idxs.len(),
    );
    let (jones_array, weight_array, flag_array) = synthesize_test_data(shape);

    crt.bench_function(
        "ms_writer::write_vis_mwalib - mwax_half_1247842824",
        |bch| {
            bch.iter(|| {
                ms_writer
                    .write_vis_mwalib(
                        jones_array.view(),
                        weight_array.view(),
                        flag_array.view(),
                        &context,
                        &mwalib_timestep_range,
                        &mwalib_coarse_chan_range,
                        &mwalib_baseline_idxs,
                    )
                    .unwrap();
            })
        },
    );
}

criterion_group!(
    name = io;
    config = Criterion::default().sample_size(NUM_ITERATIONS);
    targets =
        bench_ms_init_mwax_half_1247842824,
        bench_ms_write_mwax_part_1247842824
);

criterion_main!(io);