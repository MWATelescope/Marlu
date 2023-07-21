// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Benchmarks

#[cfg(not(all(feature = "mwalib", feature = "ms")))]
compile_error!("Benchmarks require the \"mwalib\" and \"ms\" features");

use std::{
    cmp::min,
    env,
    path::{Path, PathBuf},
};

use criterion::*;
use glob::glob;
use hifitime::Duration;
use marlu::{
    ms::MeasurementSetWriter, mwalib, ndarray::Array3, uvfits::UvfitsWriter, Complex, Jones,
    MwaObsContext, ObsContext, VisContext, VisSelection, VisWrite,
};
use mwalib::CorrelatorContext;
use tempfile::tempdir;

// ///////////// //
// IO Benchmarks //
// ///////////// //

const TIMESTEP_LIMIT: usize = 10;

fn get_test_dir() -> String {
    env::var("MARLU_TEST_DIR").unwrap_or_else(|_| String::from("/mnt/data"))
}

fn get_context_mwax_half_1247842824() -> CorrelatorContext {
    let test_dir = get_test_dir();
    let test_path = Path::new(&test_dir);
    let vis_path = test_path.join("1247842824_vis");
    let metafits_path = vis_path.join("1247842824.metafits");
    let gpufits_glob = vis_path.join("1247842824_*gpubox*_00.fits");
    let gpufits_files: Vec<PathBuf> = glob(gpufits_glob.to_str().unwrap())
        .unwrap()
        .filter_map(Result::ok)
        .collect();
    CorrelatorContext::new(metafits_path, &gpufits_files).unwrap()
}

fn bench_ms_init_mwax_half_1247842824(crt: &mut Criterion) {
    let corr_ctx = get_context_mwax_half_1247842824();

    let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
    vis_sel.timestep_range = vis_sel.timestep_range.start
        ..min(
            vis_sel.timestep_range.start + TIMESTEP_LIMIT + 1,
            vis_sel.timestep_range.end,
        );

    let obs_ctx = ObsContext::from_mwalib(&corr_ctx.metafits_context);
    let mwa_ctx = MwaObsContext::from_mwalib(&corr_ctx.metafits_context);

    let vis_ctx = VisContext::from_mwalib(
        &corr_ctx,
        &vis_sel.timestep_range,
        &vis_sel.coarse_chan_range,
        &vis_sel.baseline_idxs,
        1,
        1,
    );

    crt.bench_function(
        "MeasurementSetWriter::initialize_mwa - mwax_half_1247842824",
        |bch| {
            bch.iter(|| {
                let tmp_dir = tempdir().unwrap();
                let ms_path = tmp_dir.path().join("vis.ms");
                let ms_writer = MeasurementSetWriter::new(
                    ms_path,
                    obs_ctx.phase_centre,
                    obs_ctx.array_pos,
                    obs_ctx.ant_positions_geodetic().collect(),
                    Duration::from_total_nanoseconds(0),
                    true,
                );
                ms_writer
                    .initialize_mwa(
                        &vis_ctx,
                        &obs_ctx,
                        &mwa_ctx,
                        None,
                        &vis_sel.coarse_chan_range,
                    )
                    .unwrap();
            })
        },
    );
}

// This benchmark is a little misleading - "initialising" a uvfits takes almost
// no time, because very little is written (if anything) until visibilities are
// written. This benchmark closes the writer before any visibilities are
// written, so the amount of time to write an empty uvfits file (of the right
// size) is what is benchmarked.
fn bench_uvfits_init_mwax_half_1247842824(crt: &mut Criterion) {
    let corr_ctx = get_context_mwax_half_1247842824();

    let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
    vis_sel.timestep_range = vis_sel.timestep_range.start
        ..min(
            vis_sel.timestep_range.start + TIMESTEP_LIMIT + 1,
            vis_sel.timestep_range.end,
        );

    let obs_ctx = ObsContext::from_mwalib(&corr_ctx.metafits_context);

    let vis_ctx = VisContext::from_mwalib(
        &corr_ctx,
        &vis_sel.timestep_range,
        &vis_sel.coarse_chan_range,
        &vis_sel.baseline_idxs,
        1,
        1,
    );

    crt.bench_function(
        "UvfitsWriter::initialize_mwa - mwax_half_1247842824",
        |bch| {
            bch.iter(|| {
                let tmp_dir = tempdir().unwrap();
                let uvfits_path = tmp_dir.path().join("vis.uvfits");
                let u = UvfitsWriter::from_marlu(
                    uvfits_path,
                    &vis_ctx,
                    obs_ctx.array_pos,
                    obs_ctx.phase_centre,
                    Duration::from_total_nanoseconds(0),
                    obs_ctx.name.as_deref(),
                    obs_ctx.ant_names.clone(),
                    obs_ctx.ant_positions_geodetic().collect(),
                    true,
                    None,
                )
                .unwrap();
                u.close().unwrap();
            })
        },
    );
}

fn synthesize_test_data(
    shape: (usize, usize, usize),
) -> (Array3<Jones<f32>>, Array3<f32>, Array3<bool>) {
    let jones_array = Array3::from_shape_fn(shape, |(timestep_idx, chan_idx, baseline_idx)| {
        Jones::from([
            Complex::new(0., timestep_idx as _),
            Complex::new(0., chan_idx as _),
            Complex::new(0., baseline_idx as _),
            Complex::new(0., 1.),
        ])
    });
    let weight_array = Array3::from_shape_fn(shape, |(timestep_idx, chan_idx, baseline_idx)| {
        (timestep_idx * shape.0 * shape.1 * shape.2
            + chan_idx * shape.0 * shape.1
            + baseline_idx * shape.0) as f32
    });
    let flag_array = Array3::from_shape_fn(shape, |(timestep_idx, chan_idx, baseline_idx)| {
        (timestep_idx * shape.0 * shape.1 * shape.2
            + chan_idx * shape.0 * shape.1
            + baseline_idx * shape.0 % 2
            == 0) as _
    });
    (jones_array, weight_array, flag_array)
}

fn bench_ms_write_mwax_part_1247842824(crt: &mut Criterion) {
    let corr_ctx = get_context_mwax_half_1247842824();

    let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
    vis_sel.timestep_range = vis_sel.timestep_range.start
        ..min(
            vis_sel.timestep_range.start + TIMESTEP_LIMIT + 1,
            vis_sel.timestep_range.end,
        );

    let obs_ctx = ObsContext::from_mwalib(&corr_ctx.metafits_context);
    let mwa_ctx = MwaObsContext::from_mwalib(&corr_ctx.metafits_context);

    let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
    let shape = vis_sel.get_shape(fine_chans_per_coarse);
    let (jones_array, weight_array, _) = synthesize_test_data(shape);

    let vis_ctx = VisContext::from_mwalib(
        &corr_ctx,
        &vis_sel.timestep_range,
        &vis_sel.coarse_chan_range,
        &vis_sel.baseline_idxs,
        1,
        1,
    );

    crt.bench_function(
        &format!("MeasurementSetWriter::write_vis_marlu - mwax_half_1247842824 {shape:?}"),
        |bch| {
            bch.iter(|| {
                let tmp_dir = tempdir().unwrap();
                let ms_path = tmp_dir.path().join("vis.ms");
                let mut ms_writer = MeasurementSetWriter::new(
                    ms_path,
                    obs_ctx.phase_centre,
                    obs_ctx.array_pos,
                    obs_ctx.ant_positions_geodetic().collect(),
                    Duration::from_total_nanoseconds(0),
                    true,
                );
                ms_writer
                    .initialize_mwa(
                        &vis_ctx,
                        &obs_ctx,
                        &mwa_ctx,
                        None,
                        &vis_sel.coarse_chan_range,
                    )
                    .unwrap();
                ms_writer
                    .write_vis(jones_array.view(), weight_array.view(), &vis_ctx)
                    .unwrap();
            })
        },
    );
}

fn bench_uvfits_write_mwax_part_1247842824(crt: &mut Criterion) {
    let corr_ctx = get_context_mwax_half_1247842824();

    let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
    vis_sel.timestep_range = vis_sel.timestep_range.start
        ..min(
            vis_sel.timestep_range.start + TIMESTEP_LIMIT + 1,
            vis_sel.timestep_range.end,
        );

    let obs_ctx = ObsContext::from_mwalib(&corr_ctx.metafits_context);

    let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
    let shape = vis_sel.get_shape(fine_chans_per_coarse);
    let (jones_array, weight_array, _) = synthesize_test_data(shape);

    let vis_ctx = VisContext::from_mwalib(
        &corr_ctx,
        &vis_sel.timestep_range,
        &vis_sel.coarse_chan_range,
        &vis_sel.baseline_idxs,
        1,
        1,
    );

    crt.bench_function(
        &format!("UvfitsWriter::write_vis_marlu - mwax_half_1247842824 {shape:?}"),
        |bch| {
            bch.iter(|| {
                let tmp_dir = tempdir().unwrap();
                let uvfits_path = tmp_dir.path().join("vis.uvfits");
                let mut uvfits_writer = UvfitsWriter::from_marlu(
                    uvfits_path,
                    &vis_ctx,
                    obs_ctx.array_pos,
                    obs_ctx.phase_centre,
                    Duration::from_total_nanoseconds(0),
                    obs_ctx.name.as_deref(),
                    obs_ctx.ant_names.clone(),
                    obs_ctx.ant_positions_geodetic().collect(),
                    true,
                    None,
                )
                .unwrap();
                uvfits_writer
                    .write_vis(jones_array.view(), weight_array.view(), &vis_ctx)
                    .unwrap();
                uvfits_writer.close().unwrap();
            })
        },
    );
}

criterion_group!(
    name = io;
    config = Criterion::default().sample_size(60);
    targets =
        bench_ms_init_mwax_half_1247842824,
        bench_uvfits_init_mwax_half_1247842824,
        bench_ms_write_mwax_part_1247842824,
        bench_uvfits_write_mwax_part_1247842824,
);

criterion_main!(io);
