//! Spectral and Temporal averaging

use crate::Complex;
use itertools::izip;
use ndarray::prelude::*;
// {Array3, Array4, ArrayView3, ArrayView4, Axis};
use thiserror::Error;

use crate::Jones;

// use rubbl_casatables::CasacoreError;

#[derive(Error, Debug)]
pub enum AveragingError {
    #[error("bad array shape supplied to argument {argument} of function {function}. expected {expected}, received {received}")]
    BadArrayShape {
        argument: String,
        function: String,
        expected: String,
        received: String,
    },
    // TODO: https://github.com/pkgw/rubbl/pull/148
    // #[error("{0}")]
    // RubblError(#[from] CasacoreError)
}

/// compute the weighted geometric average of unflagged visibilities for each time, frequency and
/// pol in the chunks.
///
/// If all visibilities in the chunk are flagged, then use the geometric average withough weights.
///
/// dimensions:
/// - jones_chunk -> [at, af][pol] (each element is a complex vis array of length pol)
/// - weight_chunk -> [at, af, pol]
/// - flag_chunk shape -> [at, af, pol]
/// - avg_jones -> [pol]
/// - avg_weight_view -> [pol]
/// - avg_flag_view -> [pol]
///
/// TODO: this takes 4 flags, but nothing actually uses these, so just rip this out.
///
// average_chunk_for_pols(
// jones_chunk: ArrayView2<Jones<f32>>,
// weight_chunk: ArrayView3<f32>,
// flag_chunk: ArrayView3<bool>,
// avg_jones: &mut [Complex<f32>],
// avg_weight_view: ArrayViewMut1<f32>,
// avg_flag_view: ArrayViewMut1<bool>,
// ) {
#[macro_export]
macro_rules! average_chunk_for_pols_f64 {
    (
        // to be averaged
        $jones_chunk:expr,
        $weight_chunk:expr,
        $flag_chunk:expr,
        // to average into
        $avg_jones:expr,
        $avg_weight_view:expr,
        $avg_flag_view:expr
    ) => {
        let chunk_size = $jones_chunk.len();

        let mut weight_sum = [0_f64; 4];
        let mut jones_sum = Jones::<f64>::default();
        let mut jones_weighted_sum = Jones::<f64>::default();
        let mut all_flagged = true;

        for (jones, weight, flag) in izip!(
            $jones_chunk.axis_iter(Axis(0)),
            $weight_chunk.axis_iter(Axis(0)),
            $flag_chunk.axis_iter(Axis(0))
        ) {
            for (jones, weight, flag) in izip!(
                jones.iter(),
                weight.axis_iter(Axis(0)),
                flag.axis_iter(Axis(0))
            ) {
                let jones_c64 = Jones::<f64>::from(*jones);
                jones_sum += jones_c64;
                for (jones_elem, weight_elem, flag_elem, weighted_vis_sum, weight_sum) in izip!(
                    jones_c64.iter(),
                    weight.iter(),
                    flag.iter(),
                    jones_weighted_sum.iter_mut(),
                    weight_sum.iter_mut(),
                ) {
                    let weight_f64: f64 = *weight_elem as _;

                    if !flag_elem && *weight_elem > 0. {
                        *weighted_vis_sum += jones_elem * weight_f64;
                        *weight_sum += weight_f64;
                        all_flagged = false;
                    }
                }
            }
        }

        for (jones_weighted_sum, jones_sum, avg_weight_view, avg_jones, weight_sum) in izip!(
            jones_weighted_sum.iter(),
            jones_sum.iter(),
            $avg_weight_view.iter_mut(),
            $avg_jones.iter_mut(),
            weight_sum.iter()
        ) {
            *avg_jones = if !all_flagged {
                Complex::<f32>::new(
                    (jones_weighted_sum.re / weight_sum) as f32,
                    (jones_weighted_sum.im / weight_sum) as f32,
                )
            } else {
                Complex::<f32>::new(
                    (jones_sum.re / chunk_size as f64) as f32,
                    (jones_sum.im / chunk_size as f64) as f32,
                )
            };
            *avg_weight_view = *weight_sum as f32;
        }

        $avg_flag_view.fill(all_flagged);
    };
}

pub type VisData = (Array3<Jones<f32>>, Array4<f32>, Array4<bool>);

/// Average a section (timestep_range, coarse_chan_range) of the visibilities
/// (jones_array, weight_array, flag_array) in time or frequency (time_factor, frequency_factor).
///
/// `jones_array` - a three dimensional array of jones matrix visibilities.
///     The dimensions of the array are `[timestep][channel][baseline]`
///
/// `weight_array` - a four dimensional array of visibility weights.
///     The dimensions of the array are `[timestep][channel][baseline][pol]`
///
/// `flag_array` - a four dimensional array of visibility flags.
///     The dimensions of the array are `[timestep][channel][baseline][pol]`
///
/// `time_factor` - the factor by which to average the time axis.
///
/// `frequency_factor` - the factor by which to average the frequency axis.
///
/// # Gorey details
///
/// Averaging is done "Cotter-style". For each `time_factor` * `frequency_factor`
/// chunk of input visibilities:
/// - unflagged weights are added together
/// - if all visibilities in a chunk are flagged, then the result is the geometric
///     mean of the chunk.
/// - otherwise the visibility is the weighted mean of the unflagged visibilities.
///
/// This has been validated thoroughly against Cotter.
///
pub fn average_visibilities(
    jones_array: ArrayView3<Jones<f32>>,
    weight_array: ArrayView4<f32>,
    flag_array: ArrayView4<bool>,
    avg_time: usize,
    avg_freq: usize,
) -> Result<VisData, AveragingError> {
    let jones_dims = jones_array.dim();
    let weight_dims = weight_array.dim();
    if weight_dims != (jones_dims.0, jones_dims.1, jones_dims.2, 4) {
        return Err(AveragingError::BadArrayShape {
            argument: "weight_array".to_string(),
            function: "average_visibilities".to_string(),
            expected: format!("({}, {}, {}, 4)", jones_dims.0, jones_dims.1, jones_dims.2),
            received: format!("{:?}", weight_dims),
        });
    }
    let flag_dims = flag_array.dim();
    if flag_dims != (jones_dims.0, jones_dims.1, jones_dims.2, 4) {
        return Err(AveragingError::BadArrayShape {
            argument: "flag_array".to_string(),
            function: "average_visibilities".to_string(),
            expected: format!("({}, {}, {}, 4)", jones_dims.0, jones_dims.1, jones_dims.2),
            received: format!("{:?}", flag_dims),
        });
    }
    let averaged_dims = (
        (jones_dims.0 as f64 / avg_time as f64).ceil() as usize,
        (jones_dims.1 as f64 / avg_freq as f64).ceil() as usize,
        jones_dims.2,
    );
    let mut averaged_jones_array = Array3::<Jones<f32>>::zeros(averaged_dims);
    let mut averaged_weight_array =
        Array4::<f32>::zeros((averaged_dims.0, averaged_dims.1, averaged_dims.2, 4));
    let mut averaged_flag_array = Array4::<bool>::from_elem(
        (averaged_dims.0, averaged_dims.1, averaged_dims.2, 4),
        false,
    );

    // iterate through the time dimension of the arrays in chunks of size `time_factor`.
    for (
        jones_timestep_chunk,
        weight_timestep_chunk,
        flag_timestep_chunk,
        mut averaged_jones_timestep_view,
        mut averaged_weight_timestep_view,
        mut averaged_flag_timestep_view,
    ) in izip!(
        jones_array.axis_chunks_iter(Axis(0), avg_time),
        weight_array.axis_chunks_iter(Axis(0), avg_time),
        flag_array.axis_chunks_iter(Axis(0), avg_time),
        averaged_jones_array.outer_iter_mut(),
        averaged_weight_array.outer_iter_mut(),
        averaged_flag_array.outer_iter_mut(),
    ) {
        // iterate through the channel dimension of the arrays in chunks of size `frequency_factor`.
        for (
            jones_channel_chunk,
            weight_channel_chunk,
            flag_channel_chunk,
            mut averaged_jones_channel_view,
            mut averaged_weight_channel_view,
            mut averaged_flag_channel_view,
        ) in izip!(
            jones_timestep_chunk.axis_chunks_iter(Axis(1), avg_freq),
            weight_timestep_chunk.axis_chunks_iter(Axis(1), avg_freq),
            flag_timestep_chunk.axis_chunks_iter(Axis(1), avg_freq),
            averaged_jones_timestep_view.outer_iter_mut(),
            averaged_weight_timestep_view.outer_iter_mut(),
            averaged_flag_timestep_view.outer_iter_mut(),
        ) {
            // iterate through the baseline dimension of the arrays.
            for (
                jones_chunk,
                weight_chunk,
                flag_chunk,
                mut averaged_jones_view,
                mut averaged_weight_view,
                mut averaged_flag_view,
            ) in izip!(
                jones_channel_chunk.axis_iter(Axis(2)),
                weight_channel_chunk.axis_iter(Axis(2)),
                flag_channel_chunk.axis_iter(Axis(2)),
                averaged_jones_channel_view.outer_iter_mut(),
                averaged_weight_channel_view.outer_iter_mut(),
                averaged_flag_channel_view.outer_iter_mut(),
            ) {
                average_chunk_for_pols_f64!(
                    jones_chunk,
                    weight_chunk,
                    flag_chunk,
                    averaged_jones_view[()],
                    averaged_weight_view,
                    averaged_flag_view
                );
            }
        }
    }

    Ok((
        averaged_jones_array,
        averaged_weight_array,
        averaged_flag_array,
    ))
}

#[cfg(test)]
mod tess {
    use crate::Complex;
    use approx::assert_abs_diff_eq;
    use ndarray::{Array3, Array4};

    use super::*;

    fn synthesize_test_data(
        shape: (usize, usize, usize, usize),
    ) -> (Array3<Jones<f32>>, Array4<f32>, Array4<bool>) {
        let vis_array = Array3::from_shape_fn(
            (shape.0, shape.1, shape.2),
            |(timestep_idx, chan_idx, baseline_idx)| {
                Jones::from([
                    Complex::new(0., timestep_idx as _),
                    Complex::new(0., chan_idx as _),
                    Complex::new(0., baseline_idx as _),
                    Complex::new(0., 1.),
                ])
            },
        );
        let weight_array =
            Array4::from_shape_fn(shape, |(timestep_idx, chan_idx, baseline_idx, pol_idx)| {
                (timestep_idx * shape.0 * shape.1 * shape.2
                    + chan_idx * shape.0 * shape.1
                    + baseline_idx * shape.0
                    + pol_idx
                    + 1) as f32
            });
        let flag_array =
            Array4::from_shape_fn(shape, |(timestep_idx, chan_idx, baseline_idx, pol_idx)| {
                (timestep_idx * shape.0 * shape.1 * shape.2
                    + chan_idx * shape.0 * shape.1
                    + baseline_idx * shape.0
                    + pol_idx % 2
                    == 0) as _
            });
        (vis_array, weight_array, flag_array)
    }

    #[test]
    fn test_averaging_trivial() {
        // let n_ants = 128;
        let n_ants = 3;
        let n_timesteps = 5;
        let n_channels = 7;
        let n_pols = 4;
        let n_baselines = n_ants * (n_ants - 1) / 2;
        let shape = (n_timesteps, n_channels, n_baselines, n_pols);
        // generate synthetic test data
        let (vis_array, weight_array, flag_array) = synthesize_test_data(shape);
        let no_flags = Array4::from_elem(flag_array.dim(), false);

        // trivial case: no averaging or flags, check that nothing weird happens.
        let (averaged_vis_array, averaged_weight_array, averaged_flag_array) =
            average_visibilities(vis_array.view(), weight_array.view(), no_flags.view(), 1, 1)
                .unwrap();

        assert_eq!(averaged_vis_array.dim(), (5, 7, 3));
        assert_eq!(averaged_weight_array.dim(), (5, 7, 3, 4));
        assert_eq!(averaged_flag_array.dim(), (5, 7, 3, 4));

        assert_abs_diff_eq!(averaged_vis_array, vis_array.view());
        assert_abs_diff_eq!(averaged_weight_array, weight_array.view());
    }

    #[test]
    fn test_averaging_non_divisors() {
        // let n_ants = 128;
        let n_ants = 3;
        let n_timesteps = 5;
        let n_channels = 7;
        let n_pols = 4;
        let n_baselines = n_ants * (n_ants - 1) / 2;
        let shape = (n_timesteps, n_channels, n_baselines, n_pols);
        let (time_factor, frequency_factor) = (2, 2);
        // generate synthetic test data
        let (vis_array, weight_array, flag_array) = synthesize_test_data(shape);

        // flags, averaging factors not divisors of shape
        let (averaged_vis_array, averaged_weight_array, averaged_flag_array) =
            average_visibilities(
                vis_array.view(),
                weight_array.view(),
                flag_array.view(),
                time_factor,
                frequency_factor,
            )
            .unwrap();

        assert_eq!(averaged_vis_array.dim(), (3, 4, 3));
        assert_eq!(averaged_weight_array.dim(), (3, 4, 3, 4));
        assert_eq!(averaged_flag_array.dim(), (3, 4, 3, 4));

        #[rustfmt::skip]
        let expected_weight_0_0_0_0 =
            (if flag_array[(0, 0, 0, 0)] {0.} else {weight_array[(0, 0, 0, 0)]})
                + (if flag_array[(1, 0, 0, 0)] {0.} else {weight_array[(1, 0, 0, 0)]})
                + (if flag_array[(0, 1, 0, 0)] {0.} else {weight_array[(0, 1, 0, 0)]})
                + (if flag_array[(1, 1, 0, 0)] {0.} else {weight_array[(1, 1, 0, 0)]});
        assert_abs_diff_eq!(averaged_weight_array[(0, 0, 0, 0)], expected_weight_0_0_0_0);

        #[rustfmt::skip]
        let expected_weight_2_3_2_3 =
            if flag_array[(4, 6, 2, 3)] {0.} else {weight_array[(4, 6, 2, 3)]};
        assert_abs_diff_eq!(averaged_weight_array[(2, 3, 2, 3)], expected_weight_2_3_2_3);
    }

    // TODO: test unflagged with zero weight.
}
