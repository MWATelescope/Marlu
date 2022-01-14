//! Spectral and Temporal averaging

use crate::Complex;
use itertools::izip;
use ndarray::{ s, Array3, Array4, ArrayView3, ArrayView4, Axis };
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
/// It's up to you to initialize the avg_ views to zero;
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
macro_rules! average_chunk_for_pols {
    (
        $jones_chunk:expr,
        $weight_chunk:expr,
        $flag_chunk:expr,
        $avg_jones:expr,
        $avg_weight_view:expr,
        $avg_flag_view:expr
    ) => {
        // let chunk_size: usize = $jones_chunk.shape()[0..1].iter().product();
        let chunk_size: usize = $jones_chunk.len();

        assert_eq!($weight_chunk.dim().2, 4);
        assert_eq!($flag_chunk.dim().2, 4);
        assert_eq!($avg_weight_view.len(), 4);
        assert_eq!($avg_flag_view.len(), 4);

        for pol_idx in 0..4 {
            $avg_flag_view[pol_idx] = true;
            // $avg_jones[pol_idx] = Complex::default();
            // $avg_weight_view[pol_idx] = -0.;
            let mut weight_sum: f64 = -0.;
            let mut complex_sum_c64 = Complex::<f64>::default();
            let mut weighted_complex_sum_c64 = Complex::<f64>::default();

            for (jones_element, weight_element, flag_element) in izip!(
                $jones_chunk.iter(),
                $weight_chunk.slice(s![.., .., pol_idx]).iter(),
                $flag_chunk.slice(s![.., .., pol_idx]).iter(),
            ) {
                let pol_c64 = Complex::<f64>::new(
                    jones_element[pol_idx].re as _,
                    jones_element[pol_idx].im as _,
                );
                complex_sum_c64 += pol_c64;
                if !flag_element && weight_element > &0. {
                    weighted_complex_sum_c64 += (pol_c64 * *weight_element as f64);
                    weight_sum += *weight_element as f64;
                    $avg_flag_view[pol_idx] = false;
                }
            }

            weighted_complex_sum_c64 = if !$avg_flag_view[pol_idx] {
                weighted_complex_sum_c64 / weight_sum
            } else {
                complex_sum_c64 / (chunk_size as f64)
            };
            $avg_jones[pol_idx] = Complex::new(
                weighted_complex_sum_c64.re as _,
                weighted_complex_sum_c64.im as _,
            );
            $avg_weight_view[pol_idx] = weight_sum as f32;
        }
    };
}

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
/// # TODO:
///
/// - visibilities with weights of zero or less should be considered flagged.
/// - if one pol is flagged, should all pols be?
/// - does it only make sense to offer power of 2 averaging?
///
pub fn average_visibilities(
    jones_array: ArrayView3<Jones<f32>>,
    weight_array: ArrayView4<f32>,
    flag_array: ArrayView4<bool>,
    avg_time: usize,
    avg_freq: usize,
) -> Result<(Array3<Jones<f32>>, Array4<f32>, Array4<bool>), AveragingError> {
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
                average_chunk_for_pols!(
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
