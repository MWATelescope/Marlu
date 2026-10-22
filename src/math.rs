// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Some helper mathematics.

/// Convert a _cross-correlation_ baseline index into its constituent tile
/// indices. Baseline 0 _is not_ between tile 0 and tile 0; it is between tile 0
/// and tile 1.
// Courtesy Brian Crosse.
#[inline]
pub fn cross_correlation_baseline_to_tiles(
    total_num_tiles: usize,
    baseline: usize,
) -> (usize, usize) {
    let n = (total_num_tiles - 1) as f64;
    let bl = baseline as f64;
    let tile1 = (-0.5 * (4.0 * n * (n + 1.0) - 8.0 * bl + 1.0).sqrt() + n + 0.5).floor();
    let tile2 = bl - tile1 * (n - (tile1 + 1.0) / 2.0) + 1.0;
    (tile1 as usize, tile2 as usize)
}

/// Convert a baseline index into its constituent tile indices (where the
/// baseline indices include auto-correlations as baselines). Baseline 0 is
/// between tile 0 and tile 0.
// Courtesy Brian Crosse.
#[inline]
pub fn baseline_to_tiles(total_num_tiles: usize, baseline: usize) -> (usize, usize) {
    let n = total_num_tiles as f64;
    let bl = baseline as f64;
    let tile1 = (-0.5 * (4.0 * n * (n + 1.0) - 8.0 * bl + 1.0).sqrt() + n + 0.5).floor();
    let tile2 = bl - tile1 * (n - (tile1 + 1.0) / 2.0);
    (tile1 as usize, tile2 as usize)
}

/// Convert tile indices into a _cross-correlation_ baseline index; the inverse
/// of [`cross_correlation_baseline_to_tiles`]. `tile1` must be less than
/// `tile2`.
#[inline]
pub fn tiles_to_cross_correlation_baseline(
    total_num_tiles: usize,
    tile1: usize,
    tile2: usize,
) -> usize {
    debug_assert!(tile1 < tile2 && tile2 < total_num_tiles);
    tile1 * (2 * total_num_tiles - tile1 - 1) / 2 + tile2 - tile1 - 1
}

/// Convert tile indices into a baseline index (where the baseline indices
/// include auto-correlations as baselines); the inverse of
/// [`baseline_to_tiles`]. `tile1` must not be greater than `tile2`.
#[inline]
pub fn tiles_to_baseline(total_num_tiles: usize, tile1: usize, tile2: usize) -> usize {
    debug_assert!(tile1 <= tile2 && tile2 < total_num_tiles);
    tile1 * (2 * total_num_tiles - tile1 + 1) / 2 + tile2 - tile1
}

/// From the number of cross-correlation baselines, get the number of tiles.
// From the definition of how many baselines there are in an array of N tiles,
// this is just the solved quadratic.
#[inline]
pub fn num_tiles_from_num_cross_correlation_baselines(num_baselines: usize) -> usize {
    (((1 + 8 * num_baselines) as f64).sqrt() as usize).div_ceil(2)
}

/// From the number of baselines (which also include auto-correlations as
/// baselines), get the number of tiles.
// From the definition of how many baselines there are in an array of N tiles,
// this is just the solved quadratic.
#[inline]
pub fn num_tiles_from_num_baselines(num_baselines: usize) -> usize {
    (((1 + 8 * num_baselines) as f64).sqrt() as usize - 1) / 2
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cross_correlation_baseline_to_tiles() {
        // Let's pretend we have 128 tiles, therefore 8128 baselines. Check that
        // our function does the right thing.
        let n = 128;
        let mut bl_index = 0;
        for tile1 in 0..n {
            for tile2 in tile1 + 1..n {
                let (t1, t2) = cross_correlation_baseline_to_tiles(n, bl_index);
                assert_eq!(
                    tile1, t1,
                    "Expected tile1 = {tile1}, got {t1}. bl = {bl_index}"
                );
                assert_eq!(
                    tile2, t2,
                    "Expected tile2 = {tile2}, got {t2}. bl = {bl_index}"
                );
                bl_index += 1;
            }
        }

        // Try with a different number of tiles.
        let n = 126;
        let mut bl_index = 0;
        for tile1 in 0..n {
            for tile2 in tile1 + 1..n {
                let (t1, t2) = cross_correlation_baseline_to_tiles(n, bl_index);
                assert_eq!(
                    tile1, t1,
                    "Expected tile1 = {tile1}, got {t1}. bl = {bl_index}"
                );
                assert_eq!(
                    tile2, t2,
                    "Expected tile2 = {tile2}, got {t2}. bl = {bl_index}"
                );
                bl_index += 1;
            }
        }

        let n = 256;
        let mut bl_index = 0;
        for tile1 in 0..n {
            for tile2 in tile1 + 1..n {
                let (t1, t2) = cross_correlation_baseline_to_tiles(n, bl_index);
                assert_eq!(
                    tile1, t1,
                    "Expected tile1 = {tile1}, got {t1}. bl = {bl_index}"
                );
                assert_eq!(
                    tile2, t2,
                    "Expected tile2 = {tile2}, got {t2}. bl = {bl_index}"
                );
                bl_index += 1;
            }
        }
    }

    #[test]
    fn test_tiles_to_baselines() {
        for n in [1, 2, 3, 126, 128, 256] {
            let mut cross = 0;
            let mut all = 0;
            for tile1 in 0..n {
                for tile2 in tile1..n {
                    assert_eq!(tiles_to_baseline(n, tile1, tile2), all);
                    assert_eq!(baseline_to_tiles(n, all), (tile1, tile2));
                    all += 1;
                    if tile1 != tile2 {
                        assert_eq!(tiles_to_cross_correlation_baseline(n, tile1, tile2), cross);
                        assert_eq!(
                            cross_correlation_baseline_to_tiles(n, cross),
                            (tile1, tile2)
                        );
                        cross += 1;
                    }
                }
            }
            assert_eq!(all, n * (n + 1) / 2);
            assert_eq!(cross, n * (n - 1) / 2);
        }
    }

    #[test]
    fn test_baseline_to_tiles() {
        // Let's pretend we have 128 tiles, therefore 8256 baselines. Check that
        // our function does the right thing.
        let n = 128;
        let mut bl_index = 0;
        for tile1 in 0..n {
            for tile2 in tile1..n {
                let (t1, t2) = baseline_to_tiles(n, bl_index);
                assert_eq!(
                    tile1, t1,
                    "Expected tile1 = {tile1}, got {t1}. bl = {bl_index}"
                );
                assert_eq!(
                    tile2, t2,
                    "Expected tile2 = {tile2}, got {t2}. bl = {bl_index}"
                );
                bl_index += 1;
            }
        }

        // Try with a different number of tiles.
        let n = 126;
        let mut bl_index = 0;
        for tile1 in 0..n {
            for tile2 in tile1..n {
                let (t1, t2) = baseline_to_tiles(n, bl_index);
                assert_eq!(
                    tile1, t1,
                    "Expected tile1 = {tile1}, got {t1}. bl = {bl_index}"
                );
                assert_eq!(
                    tile2, t2,
                    "Expected tile2 = {tile2}, got {t2}. bl = {bl_index}"
                );
                bl_index += 1;
            }
        }

        let n = 256;
        let mut bl_index = 0;
        for tile1 in 0..n {
            for tile2 in tile1..n {
                let (t1, t2) = baseline_to_tiles(n, bl_index);
                assert_eq!(
                    tile1, t1,
                    "Expected tile1 = {tile1}, got {t1}. bl = {bl_index}"
                );
                assert_eq!(
                    tile2, t2,
                    "Expected tile2 = {tile2}, got {t2}. bl = {bl_index}"
                );
                bl_index += 1;
            }
        }
    }

    #[test]
    fn test_num_tiles_from_num_cross_correlation_baselines() {
        assert_eq!(num_tiles_from_num_cross_correlation_baselines(8128), 128);
        assert_eq!(num_tiles_from_num_cross_correlation_baselines(8001), 127);
        assert_eq!(num_tiles_from_num_cross_correlation_baselines(15), 6);
    }

    #[test]
    fn test_num_tiles_from_num_baselines() {
        assert_eq!(num_tiles_from_num_baselines(8256), 128);
        assert_eq!(num_tiles_from_num_baselines(8128), 127);
        assert_eq!(num_tiles_from_num_baselines(21), 6);
    }
}
