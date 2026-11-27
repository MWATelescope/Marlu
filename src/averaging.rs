// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Spectral and Temporal averaging

use crate::Complex;
use itertools::izip;
use ndarray::prelude::*;
use thiserror::Error;

use crate::Jones;

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
/// - `jones_chunk` -> [at, af][pol] (each element is a complex vis array of length pol)
/// - `weight_chunk` -> [at, af, pol]
/// - `flags_chunk` -> [at, af, pol]
/// - `avg_jones` -> [pol]
/// - `avg_weight_view` -> [pol]
/// - `avg_flag_view` -> [pol]
#[macro_export]
macro_rules! average_chunk_for_pols_f64 {
    (
        // to be averaged
        $jones_chunk:expr,
        $weights_chunk:expr,
        $flags_chunk:expr,
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

        for (jones_chunk, weight_chunk, flag) in izip!(
            $jones_chunk.axis_iter(Axis(0)),
            $weights_chunk.axis_iter(Axis(0)),
            $flags_chunk.axis_iter(Axis(0))
        ) {
            for (jones, weight, flag) in izip!(
                jones_chunk.iter(),
                weight_chunk.axis_iter(Axis(0)),
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

                    if !flag_elem && *weight_elem >= 0. {
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

/// compute the weighted geometric average of unflagged visibilities for each time, and frequency
/// in the chunks. Only weights are required, instead of flags
///
/// # Weight convention
///
/// The sign of a weight is the flag: a negative weight means the visibility is
/// flagged, and its magnitude is the weight the visibility would have if it
/// were not. A weight of zero (either sign), or one that is not finite, means
/// there is no data at all, for example a timestep that is missing from the
/// raw files and was filled with zeros. Such a visibility counts as flagged,
/// is left out of every weighted average, and is only used in the plain
/// average when nothing else is available. See [`is_weight_unflagged`] and
/// [`normalise_weight`].
///
/// The average of a chunk is
///
/// - if any visibility is unflagged: the average of the unflagged
///   visibilities weighted by their weights. The result is unflagged, and its
///   weight is the sum of those weights. Flagged visibilities are ignored.
/// - otherwise, if any flagged visibility has a weight: the average of the
///   flagged visibilities weighted by the magnitudes of their weights. The
///   result is flagged, and its weight is minus the sum of those magnitudes,
///   so that it is written the same way a flagged visibility is when there is
///   no averaging at all. Visibilities with no weight don't dilute the others
///   (Birli issue #162).
/// - otherwise (every weight is zero, or not finite): the plain average of
///   the visibilities. The result is flagged, with weight `-0.0`.
///
/// The flag is also returned on its own, so a caller doesn't have to look at
/// the sign of a zero.
///
/// dimensions:
/// - `jones_chunk` -> [at, af][pol] (each element is a complex vis array of length pol)
/// - `weight_chunk` -> [at, af]
/// - `flag_chunk` -> [at, af]
/// - `avg_jones` -> [pol]
/// - `avg_weight` -> (scalar)
/// - `avg_flag` -> (scalar)
#[macro_export]
macro_rules! average_chunk_f64 {
    (
        // to be averaged
        $jones_chunk:expr,
        $weights_chunk:expr,
        // to average into
        $avg_jones:expr,
        $avg_weight:expr,
        $avg_flag:expr
    ) => {
        let (avg_jones, avg_weight, avg_flag) =
            $crate::averaging::average_chunk($jones_chunk.view(), $weights_chunk.view());
        for (from, to) in avg_jones.iter().zip($avg_jones.iter_mut()) {
            *to = *from;
        }
        $avg_weight = avg_weight;
        $avg_flag = avg_flag;
        // not every caller reads the flag
        let _ = $avg_flag;
    };
}

/// Is this the weight of an unflagged visibility? Only a finite, strictly
/// positive weight is. A negative weight is a flag, and zero (either sign) or
/// a weight that is not a number or is infinite means there is no data.
///
/// See the weight convention in [`average_chunk_f64`].
#[must_use]
pub fn is_weight_unflagged(weight: f32) -> bool {
    weight.is_finite() && weight > 0.
}

/// Write a weight the canonical way: unflagged weights are returned as they
/// are, and all others are returned as flagged, that is with a negative sign,
/// and keep their magnitude. A weight of zero becomes `-0.0` and a weight that
/// is not finite also becomes `-0.0`, since it can't be used as a weight.
///
/// Both the sign of [`normalise_weight`]`(w)` and [`is_weight_unflagged`]`(w)`
/// give the flag of a visibility with weight `w`, and they always agree.
#[must_use]
pub fn normalise_weight(weight: f32) -> f32 {
    if is_weight_unflagged(weight) {
        weight
    } else if weight.is_finite() {
        -weight.abs()
    } else {
        -0.
    }
}

/// The function behind [`average_chunk_f64`]: the average visibility, weight
/// and flag of a `[time, freq]` chunk of visibilities and weights. The weight
/// is negative if the average is flagged. See the weight convention in
/// [`average_chunk_f64`].
///
/// # Panics
///
/// Panics if the two arrays have different shapes.
pub fn average_chunk(
    jones_chunk: ArrayView2<Jones<f32>>,
    weights_chunk: ArrayView2<f32>,
) -> (Jones<f32>, f32, bool) {
    assert_eq!(
        jones_chunk.shape(),
        weights_chunk.shape(),
        "jones and weight arrays must have the same shape"
    );

    let mut weight_sum = 0_f64;
    let mut jones_sum = Jones::<f64>::default();
    let mut jones_weighted_sum = Jones::<f64>::default();
    // for the flagged samples: the sum of the weights' magnitudes, and
    // the visibilities weighted by them.
    let mut flagged_weight_sum = 0_f64;
    let mut jones_flagged_sum = Jones::<f64>::default();

    for (jones, weight) in izip!(jones_chunk.iter(), weights_chunk.iter()) {
        let jones_c64 = Jones::<f64>::from(*jones);
        let weight_abs = (*weight as f64).abs();
        jones_sum += jones_c64;
        // a weight that is zero or not finite is used by no weighted average
        if is_weight_unflagged(*weight) {
            weight_sum += weight_abs;
            jones_weighted_sum += jones_c64 * weight_abs;
        } else if weight.is_finite() && weight_abs > 0. {
            flagged_weight_sum += weight_abs;
            jones_flagged_sum += jones_c64 * weight_abs;
        }
    }

    let avg_flag = weight_sum <= 0.;
    let (avg_jones, avg_weight) = if !avg_flag {
        (jones_weighted_sum / weight_sum, weight_sum as f32)
    } else if flagged_weight_sum > 0. {
        (
            jones_flagged_sum / flagged_weight_sum,
            -(flagged_weight_sum as f32),
        )
    } else {
        (jones_sum / jones_chunk.len() as f64, -0.)
    };

    (Jones::<f32>::from(avg_jones), avg_weight, avg_flag)
}

pub type VisData344 = (Array3<Jones<f32>>, Array4<f32>, Array4<bool>);
pub type VisData33 = (Array3<Jones<f32>>, Array3<f32>);

/// Average a section (`timestep_range`, `coarse_chan_range`) of the visibilities
/// (`jones_array`, `weight_array`, `flag_array`) in time or frequency (`time_factor`, `frequency_factor`).
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
///   mean of the chunk.
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
) -> Result<VisData344, AveragingError> {
    let jones_dims = jones_array.dim();
    let weight_dims = weight_array.dim();
    if weight_dims != (jones_dims.0, jones_dims.1, jones_dims.2, 4) {
        return Err(AveragingError::BadArrayShape {
            argument: "weight_array".to_string(),
            function: "average_visibilities".to_string(),
            expected: format!("({}, {}, {}, 4)", jones_dims.0, jones_dims.1, jones_dims.2),
            received: format!("{weight_dims:?}"),
        });
    }
    let flag_dims = flag_array.dim();
    if flag_dims != (jones_dims.0, jones_dims.1, jones_dims.2, 4) {
        return Err(AveragingError::BadArrayShape {
            argument: "flag_array".to_string(),
            function: "average_visibilities".to_string(),
            expected: format!("({}, {}, {}, 4)", jones_dims.0, jones_dims.1, jones_dims.2),
            received: format!("{flag_dims:?}"),
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

/// Visibilities and weights, and what averaging them should give, shared by the
/// tests of the writers that average, so that they all follow the same weight
/// convention (see [`average_chunk_f64`]) with and without averaging.
#[cfg(test)]
pub(crate) mod test_scenarios {
    use super::{is_weight_unflagged, normalise_weight};
    use crate::Jones;
    use ndarray::Array3;

    /// A 2 timestep by 2 channel chunk of visibilities (each of which is
    /// `Jones::identity() * vis`), their weights, and the average of the
    /// chunk. The sample at timestep `t` and channel `c` is `[t][c % 2]`.
    pub(crate) struct Scenario {
        pub name: &'static str,
        pub vis: [[f32; 2]; 2],
        pub weights: [[f32; 2]; 2],
        /// The average of the visibilities, and the weight of the average,
        /// signed so that negative (including `-0.0`) is flagged.
        pub avg_vis: f32,
        pub avg_weight: f32,
    }

    impl Scenario {
        /// Is the average flagged?
        pub(crate) fn avg_flag(&self) -> bool {
            self.avg_weight.is_sign_negative()
        }
    }

    /// Every case: one scenario per baseline.
    pub(crate) fn scenarios() -> Vec<Scenario> {
        let nan = f32::NAN;
        let inf = f32::INFINITY;
        vec![
            Scenario {
                name: "unflagged",
                vis: [[42.; 2]; 2],
                weights: [[8.; 2]; 2],
                avg_vis: 42.,
                avg_weight: 32.,
            },
            Scenario {
                name: "flagged",
                vis: [[42.; 2]; 2],
                weights: [[-8.; 2]; 2],
                avg_vis: 42.,
                avg_weight: -32.,
            },
            // Birli issue 162: the second timestep is missing from the raw files
            Scenario {
                name: "flagged, then missing",
                vis: [[42., 42.], [0., 0.]],
                weights: [[-8., -8.], [0., 0.]],
                avg_vis: 42.,
                avg_weight: -16.,
            },
            Scenario {
                name: "flagged, then missing with negative zero",
                vis: [[42., 42.], [0., 0.]],
                weights: [[-8., -8.], [-0., -0.]],
                avg_vis: 42.,
                avg_weight: -16.,
            },
            Scenario {
                name: "mixed",
                vis: [[10., 99.], [20., 99.]],
                weights: [[8., -8.], [8., -8.]],
                avg_vis: 15.,
                avg_weight: 16.,
            },
            Scenario {
                name: "unflagged with zero weights",
                vis: [[10., 1000.], [20., 1000.]],
                weights: [[8., 0.], [8., -0.]],
                avg_vis: 15.,
                avg_weight: 16.,
            },
            Scenario {
                name: "all zero weights",
                vis: [[2., 4.], [6., 8.]],
                weights: [[0.; 2]; 2],
                avg_vis: 5.,
                avg_weight: -0.,
            },
            Scenario {
                name: "weights that are not finite",
                vis: [[1., 2.], [3., 4.]],
                weights: [[nan, inf], [-inf, -8.]],
                avg_vis: 4.,
                avg_weight: -8.,
            },
            Scenario {
                name: "flagged with different weights",
                vis: [[2., 6.], [2., 6.]],
                weights: [[-1., -3.], [-1., -3.]],
                avg_vis: 5.,
                avg_weight: -8.,
            },
        ]
    }

    /// Visibilities and weights of shape `[timestep][channel][baseline]`,
    /// where baseline `b` is `scenarios[b]`.
    pub(crate) fn arrays(
        scenarios: &[Scenario],
        num_timesteps: usize,
        num_chans: usize,
    ) -> (Array3<Jones<f32>>, Array3<f32>) {
        let shape = (num_timesteps, num_chans, scenarios.len());
        let vis = Array3::from_shape_fn(shape, |(t, c, b)| {
            Jones::<f32>::identity() * scenarios[b].vis[t % 2][c % 2]
        });
        let weights = Array3::from_shape_fn(shape, |(t, c, b)| scenarios[b].weights[t % 2][c % 2]);
        (vis, weights)
    }

    /// What a sample should be written as when there is no averaging: a flag,
    /// and the weight the sign of which is that flag.
    pub(crate) fn unaveraged(weight: f32) -> (bool, f32) {
        let normal = normalise_weight(weight);
        assert_eq!(normal.is_sign_negative(), !is_weight_unflagged(weight));
        (normal.is_sign_negative(), normal)
    }
}

#[cfg(test)]
mod tess {
    use crate::Complex;
    use approx::assert_abs_diff_eq;
    use ndarray::prelude::*;

    use super::{
        average_chunk, average_visibilities, is_weight_unflagged, normalise_weight, Jones,
    };

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

    #[test]
    /// birli issue 162 <https://github.com/MWATelescope/Birli/issues/162>
    fn test_averaging_birli_162() {
        let n_timesteps = 4;
        let n_channels = 64;
        let n_pols = 4;
        let n_baselines = 1;
        let shape = (n_timesteps, n_channels, n_baselines, n_pols);
        let time_factor = 2;
        let frequency_factor = 1;

        let jones_id: Jones<f32> = Jones::identity();
        let mut vis_array = Array3::from_elem((n_timesteps, n_channels, n_baselines), jones_id);
        let mut weight_array = Array4::from_elem(shape, 1_f32);
        let mut flag_array = Array4::from_elem(shape, false);

        // simulate missing timestep idx 3, second half of channels
        let half_chan = n_channels / 2;
        vis_array
            .slice_mut(s![3, half_chan.., 0])
            .fill(jones_id * -99999.0);
        weight_array
            .slice_mut(s![3, half_chan.., 0, ..])
            .fill(-99999.0);
        flag_array.slice_mut(s![3, half_chan.., 0, ..]).fill(true);

        let (averaged_vis_array, averaged_weight_array, averaged_flag_array) =
            average_visibilities(
                vis_array.view(),
                weight_array.view(),
                flag_array.view(),
                time_factor,
                frequency_factor,
            )
            .unwrap();

        // despite the missing value, it should still average to identity.
        averaged_vis_array.for_each(|&v| assert_abs_diff_eq!(v, jones_id));
        averaged_weight_array
            .slice(s![.., ..half_chan, 0, ..])
            .for_each(|&v| assert_abs_diff_eq!(v, 2.0));
        averaged_weight_array
            .slice(s![0, half_chan.., 0, ..])
            .for_each(|&v| assert_abs_diff_eq!(v, 2.0));
        averaged_weight_array
            .slice(s![1, half_chan.., 0, ..])
            .for_each(|&v| assert_abs_diff_eq!(v, 1.0));
        averaged_flag_array.for_each(|&v| assert!(!v));
    }

    /// Weights that are used as weights, and weights that are not, with the
    /// flag each gives.
    fn weights_and_flags() -> Vec<(f32, bool)> {
        vec![
            (8., false),
            (f32::MIN_POSITIVE, false),
            (-8., true),
            (-f32::MIN_POSITIVE, true),
            (0., true),
            (-0., true),
            (f32::NAN, true),
            (-f32::NAN, true),
            (f32::INFINITY, true),
            (f32::NEG_INFINITY, true),
        ]
    }

    /// Only a finite weight above zero is the weight of an unflagged
    /// visibility.
    #[test]
    fn test_is_weight_unflagged() {
        for (weight, flagged) in weights_and_flags() {
            assert_eq!(is_weight_unflagged(weight), !flagged, "weight {weight}");
        }
    }

    /// Normalising a weight makes its sign the flag, keeps the magnitude of
    /// weights that have one, and does nothing to a weight already normal.
    #[test]
    fn test_normalise_weight() {
        for (weight, flagged) in weights_and_flags() {
            let normal = normalise_weight(weight);
            assert_eq!(normal.is_sign_negative(), flagged, "weight {weight}");
            // idempotent, with the sign of a zero kept
            assert_eq!(
                normalise_weight(normal).to_bits(),
                normal.to_bits(),
                "weight {weight}"
            );
            if weight.is_finite() {
                assert_eq!(normal.abs(), weight.abs());
            } else {
                assert_eq!(normal.to_bits(), (-0_f32).to_bits());
            }
        }
        // a zero is flagged whichever way it is signed
        assert_eq!(normalise_weight(0.).to_bits(), (-0_f32).to_bits());
        assert_eq!(normalise_weight(-0.).to_bits(), (-0_f32).to_bits());
        // flagged and unflagged weights are not changed
        assert_eq!(normalise_weight(8.), 8.);
        assert_eq!(normalise_weight(-8.), -8.);
    }

    /// Average a [time, freq] chunk with [`average_chunk_f64`].
    fn average_chunk_macro(
        vis: &Array2<Jones<f32>>,
        weights: &Array2<f32>,
    ) -> (Jones<f32>, f32, bool) {
        let mut avg_jones = Jones::<f32>::default();
        let avg_weight: f32;
        let avg_flag: bool;
        average_chunk_f64!(vis, weights, avg_jones, avg_weight, avg_flag);
        (avg_jones, avg_weight, avg_flag)
    }

    /// Average a chunk both ways, check they agree, and that the sign of the
    /// weight is the flag.
    fn average_chunk_checked(
        vis: &Array2<Jones<f32>>,
        weights: &Array2<f32>,
    ) -> (Jones<f32>, f32, bool) {
        let direct = average_chunk(vis.view(), weights.view());
        let (jones, weight, flag) = average_chunk_macro(vis, weights);
        assert_abs_diff_eq!(jones, direct.0);
        assert_eq!(weight.to_bits(), direct.1.to_bits());
        assert_eq!(flag, direct.2);
        assert_eq!(weight.is_sign_negative(), flag, "weight {weight}");
        direct
    }

    /// Unflagged visibilities are averaged by weight, and flagged ones are
    /// left out.
    #[test]
    fn test_average_chunk_ignores_flagged() {
        let id = Jones::<f32>::identity();
        let vis = array![[id * 2., id * 4.], [id * 100., id * 6.]];
        let weights = array![[1., 3.], [-8., 0.]];
        let (avg_jones, avg_weight, avg_flag) = average_chunk_checked(&vis, &weights);
        assert_abs_diff_eq!(avg_jones, id * 3.5);
        assert_abs_diff_eq!(avg_weight, 4.);
        assert!(!avg_flag);
    }

    /// An unflagged visibility with zero weight is no data, so it is left out
    /// of the average (resolves an old TODO), whichever way the zero is signed.
    #[test]
    fn test_average_chunk_ignores_zero_weight_among_unflagged() {
        let id = Jones::<f32>::identity();
        let vis = array![[id * 2., id * 1000., id * 1000.]];
        for zero in [0., -0.] {
            let weights = array![[5., zero, zero]];
            let (avg_jones, avg_weight, avg_flag) = average_chunk_checked(&vis, &weights);
            assert_abs_diff_eq!(avg_jones, id * 2.);
            assert_abs_diff_eq!(avg_weight, 5.);
            assert!(!avg_flag);
        }
    }

    /// Weights that are not finite are not weights, in an unflagged average
    /// as much as in a flagged one.
    #[test]
    fn test_average_chunk_ignores_non_finite_weights() {
        let id = Jones::<f32>::identity();
        let vis = array![[id * 2., id * 1000., id * 1000., id * 1000.]];
        for bad in [f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
            let weights = array![[5., bad, bad, bad]];
            let (avg_jones, avg_weight, avg_flag) = average_chunk_checked(&vis, &weights);
            assert_abs_diff_eq!(avg_jones, id * 2.);
            assert_abs_diff_eq!(avg_weight, 5.);
            assert!(!avg_flag);

            // flagged, so the unflagged one is not a weight either
            let weights = array![[-5., bad, bad, bad]];
            let (avg_jones, avg_weight, avg_flag) = average_chunk_checked(&vis, &weights);
            assert_abs_diff_eq!(avg_jones, id * 2.);
            assert_abs_diff_eq!(avg_weight, -5.);
            assert!(avg_flag);
        }
    }

    /// birli issue 162 <https://github.com/MWATelescope/Birli/issues/162>:
    /// a flagged visibility averaged with a missing one (zeros, with no
    /// weight) keeps its amplitude.
    #[test]
    fn test_average_chunk_all_flagged_skips_zero_weight() {
        let id = Jones::<f32>::identity();
        let vis = array![[id * 42.], [id * 0.]];
        for missing_weight in [0., -0.] {
            let weights = array![[-8.], [missing_weight]];
            let (avg_jones, avg_weight, avg_flag) = average_chunk_checked(&vis, &weights);
            assert_abs_diff_eq!(avg_jones, id * 42.);
            // the weight is the flagged one: what the flagged sample would
            // have had
            assert_abs_diff_eq!(avg_weight, -8.);
            assert!(avg_flag);
        }
    }

    /// When every visibility is flagged, they are averaged by the weights
    /// they would have had, and the weights add up, as they do when they are
    /// not flagged.
    #[test]
    fn test_average_chunk_all_flagged_weighted() {
        let id = Jones::<f32>::identity();
        let vis = array![[id * 2., id * 6.]];
        let (avg_jones, avg_weight, avg_flag) = average_chunk_checked(&vis, &array![[-8., -8.]]);
        assert_abs_diff_eq!(avg_jones, id * 4.);
        assert_abs_diff_eq!(avg_weight, -16.);
        assert!(avg_flag);

        let (avg_jones, avg_weight, avg_flag) = average_chunk_checked(&vis, &array![[-1., -3.]]);
        assert_abs_diff_eq!(avg_jones, id * 5.);
        assert_abs_diff_eq!(avg_weight, -4.);
        assert!(avg_flag);

        // a weight that is not a number is not a weight
        let (avg_jones, avg_weight, avg_flag) =
            average_chunk_checked(&vis, &array![[-1., f32::NAN]]);
        assert_abs_diff_eq!(avg_jones, id * 2.);
        assert_abs_diff_eq!(avg_weight, -1.);
        assert!(avg_flag);
    }

    /// With no weights at all, the result is the plain average, and flagged,
    /// whichever way the zeros are signed, or even if they are not numbers.
    #[test]
    fn test_average_chunk_all_zero_weight() {
        let id = Jones::<f32>::identity();
        let vis = array![[id * 2., id * 6.]];
        for weights in [
            array![[0., 0.]],
            array![[-0., -0.]],
            array![[0., -0.]],
            array![[f32::NAN, f32::INFINITY]],
            array![[f32::NEG_INFINITY, 0.]],
        ] {
            let (avg_jones, avg_weight, avg_flag) = average_chunk_checked(&vis, &weights);
            assert_abs_diff_eq!(avg_jones, id * 4.);
            assert_eq!(avg_weight.to_bits(), (-0_f32).to_bits());
            assert!(avg_flag);
        }
    }

    /// A chunk with one visibility is that visibility, and the weight and
    /// flag are what the writers write when there is no averaging. This keeps
    /// the output the same whether or not it is averaged.
    #[test]
    fn test_average_chunk_of_one_is_not_averaging() {
        let jones = Jones::<f32>::identity() * 42.;
        for (weight, flagged) in weights_and_flags() {
            let (avg_jones, avg_weight, avg_flag) =
                average_chunk_checked(&array![[jones]], &array![[weight]]);
            assert_abs_diff_eq!(avg_jones, jones);
            assert_eq!(avg_weight.to_bits(), normalise_weight(weight).to_bits());
            assert_eq!(avg_flag, flagged, "weight {weight}");
        }
    }

    /// The flag of the average is whether any weight is unflagged, and it is
    /// the sign of the weight, for every combination of weights.
    #[test]
    fn test_average_chunk_flag_combinations() {
        let id = Jones::<f32>::identity();
        let vis = array![[id * 1., id * 2.]];
        for (w1, f1) in weights_and_flags() {
            for (w2, f2) in weights_and_flags() {
                let (_, avg_weight, avg_flag) = average_chunk_checked(&vis, &array![[w1, w2]]);
                // flagged only if both are
                assert_eq!(avg_flag, f1 && f2, "weights {w1} and {w2}");
                // not flagged: the weights of the unflagged ones add. flagged:
                // minus the magnitudes of the finite ones.
                let expected: f32 = [(w1, f1), (w2, f2)]
                    .iter()
                    .filter(|(_, f)| *f == avg_flag)
                    .filter(|(w, _)| w.is_finite())
                    .map(|(w, _)| w.abs())
                    .sum();
                assert_abs_diff_eq!(
                    avg_weight.abs(),
                    expected,
                    epsilon = f32::EPSILON * expected
                );
            }
        }
    }

    /// Averaging time and frequency in one go gives the same answer as
    /// averaging each in turn, when nothing is flagged, and when everything
    /// is.
    #[test]
    fn test_average_chunk_2d_chunk_layout() {
        let id = Jones::<f32>::identity();
        let vis = array![[id * 1., id * 2.], [id * 3., id * 4.]];
        for sign in [1., -1.] {
            let weights = array![[1., 2.], [3., 4.]] * sign;
            let (avg_jones, avg_weight, avg_flag) = average_chunk_checked(&vis, &weights);
            assert_abs_diff_eq!(avg_jones, id * 3.);
            assert_abs_diff_eq!(avg_weight, 10. * sign);
            assert_eq!(avg_flag, sign < 0.);
        }
    }

    /// averaging a whole array flags the output samples that have no
    /// unflagged visibility, and not the others, just like for a single chunk.
    #[test]
    fn test_average_visibilities_flagged_baseline_matches_trivial() {
        let id = Jones::<f32>::identity();
        // [time, freq, baseline, pol]: 2 timesteps, 2 channels, 2 baselines.
        let vis = Array3::from_shape_fn((2, 2, 2), |(t, f, b)| id * (1 + t + f + b) as f32);
        // baseline 1 is all flagged, with weights of magnitude 8
        let weights = Array4::from_shape_fn(
            (2, 2, 2, 4),
            |(_, _, b, _)| {
                if b == 1 {
                    -8_f32
                } else {
                    8_f32
                }
            },
        );
        let flags = weights.map(|w| *w < 0.);

        let (_, avg_weights, avg_flags) =
            average_visibilities(vis.view(), weights.view(), flags.view(), 2, 2).unwrap();
        // 4 samples averaged: baseline 0 is unflagged and baseline 1 is flagged
        for (b, expected_flag) in [(0, false), (1, true)] {
            assert_eq!(avg_flags[(0, 0, b, 0)], expected_flag);
        }
        assert_abs_diff_eq!(avg_weights[(0, 0, 0, 0)], 32.);
        // flagged, so no weight is left
        assert_abs_diff_eq!(avg_weights[(0, 0, 1, 0)], 0.);
    }

    /// The scenarios the writers are tested with are what [`average_chunk`]
    /// gives.
    #[test]
    fn test_scenarios_are_averages() {
        for scenario in super::test_scenarios::scenarios() {
            let vis = Array2::from_shape_fn((2, 2), |(t, c)| {
                Jones::<f32>::identity() * scenario.vis[t][c]
            });
            let weights = Array2::from_shape_fn((2, 2), |(t, c)| scenario.weights[t][c]);
            let (avg_jones, avg_weight, avg_flag) = average_chunk_checked(&vis, &weights);
            assert_abs_diff_eq!(
                avg_jones,
                Jones::<f32>::identity() * scenario.avg_vis,
                epsilon = 1e-5
            );
            assert_eq!(
                avg_weight.to_bits(),
                scenario.avg_weight.to_bits(),
                "{}",
                scenario.name
            );
            assert_eq!(avg_flag, scenario.avg_flag(), "{}", scenario.name);
        }
    }
}
