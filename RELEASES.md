<!-- markdownlint-disable=MD025 -->

# Version 0.12.0 (2024-08-14)

- fix issues compiling on arm64:
  - update rubbl 0.8.0 (which uses casacore v3.5.0)
  - update ndarray 0.16.0

# Version 0.11.0 (2024-05-24)

- use built 0.7, to avoid an issue in dependent crates where
  built can't find Cargo lock.
- additional error enums
- update mwalib 1.3.3

# Version 0.10.1 (2023-08-11)

- When writing out measurement sets, a weight of -0.0 is now considered a flag,
  rather than only values < 0.0 (-0.0 is not less than 0.0).

# Version 0.10.0 (2023-07-21)

- Allow vis writers to not precess their UVWs
- Improve uvfits time precision
  - A second DATE group param is now used
  - INTTIM is also used, if the time resolution was supplied
- Remove progress bars
- Remove mwalib-reading functions

# Version 0.9.2 (2023-07-18)

- update modtime when writing ms

# Version 0.9.1 (2023-02-28)

- `RADec::weighted_average` was incorrect and has now been fixed.

# Version 0.9.0 (2023-02-17)

- Change measurement sets from conditionally writing UT1 or UTC reference frames
  to always writing UTC frames. DUT1 is reported as the UT1UTC key.
- Fix a heap of clippy lints
- Remove CUDA convenience code for Rust callers
- Use mwalib v0.16.0 and fitsio v0.20.0
- Add the cargo-semver-checks action to CI
- Fix a bug in `RADec::weighted_average`
- Speed up XYZ related code
- Use the pure-Rust erfa crate rather than erfa-sys
- Rename coordinate "new" methods to "from", e.g. `RADec::new` is now
  `RADec::from_radians`
- Rename `LatLngHeight::new_mwa` to `LatLngHeight::mwa`
- Use the newest version of hifitime

# Version 0.8.0 (2022-08-22)

- Bump dependency versions.
- `cargo` feature changes:
  - `io` no longer exists
  - `cfitsio` now exists. uvfits writing is possible with just `cfitsio`, and
    `mwalib` depends on `cfitsio`.
  - `ms` now exists. Measurement Set writing is only possible with `ms`.
- Support DUT1 usage:
  - The precession API has changed
  - `UvfitsWriter` reports the DUT1 with `UT1UTC`
  - `MeasurementSetWriter` changes its time frame from `UTC` to `UT1`, iff the
    supplied DUT1 is non zero.
- IO code changes:
  - Rename `VisReadable` to `VisRead`
  - Rename `VisWritable` to `VisWrite`
  - Remove `write_vis_mwalib`
  - Rename `write_vis_marlu` to `write_vis`
  - Add a `finalise` method to `VisWrite`
  - The `UvfitsWriter` API is slightly different
- Remove a bunch of needless `clone`s from the code. This may improve
  performance.
- Add an optional `approx` feature that exposes trait implementations like
  `approx::AbsDiffEq` on each of the coordinate types (e.g. `UVW`).
- Make `Jones` `#[repr(transparent)]`

# Version 0.7.2 (2022-08-04)

- Expose Marlu version in `built_info`
- fix a bug that caused vis_ctx.timeseries to give an additional timestep.

# Version 0.7.1 (2022-08-03)

- Re-export `LmnRime`.
- Add `to_earth` and `to_earth_wgs84` functions on `XyzGeocentric`.

# Version 0.7.0 (2022-06-24)

- ⚡ @cjordan 's lightning fast uvfits optimization: using raw cfitsio instead of fitsio_sys
- use rust 1.60
- Use erfa-sys 2.0
- use ndarray 0.15.4 (instead of a range of versions)
- use mwalib 0.15.0:
  - cable lengths applied
  - expose DUT1 from metafits
- use mwalib antennas instead of rfinputs
- better error messages when creating measurement sets in paths that either don't
  exist, or are not a directory.
- api changes:
  - io:
    - uvfits `obs_name` from `Option<String>` to `Option<&str>`.
    - `history` metadata in ms and uvfits
  - Jones: convenience methods for array access
  - constants: ecpose `FREQ_WEIGHT_FACTOR`, `TIME_WEIGHT_FACTOR`
  - context: impl `Clone` for `ObsContext`
  - pos/lmn: add `LmnRime` and `LMN::prepare_for_rime`

# Version 0.6.1 (2022-03-24)

- impl Clone for VisContext
- impl PartialEq for LatLngHeight

# Version 0.6.0 (2022-03-24)

- implement VisContext, ObsContext, MwaObsContext
- migrate io::uvfits from Birli
- better error handling in io
- bake flags into weights

# Version 0.5.0 (2022-02-11)

- use mwalib v0.13.0
- kill ::time with latest hifitime
- bump min rust version from 1.56 to 1.7
- set minimum dependency versions for all deps
- Jones::nan() is more... NaNny

# Version 0.4.0 (2022-01-27)

- MeasurementSetWriter keeps track of the current row in the main table, so that rows can be written in chunks.
- MeasurementSetWriter::initialize_from_mwalib now takes the baseline_idxs array so that it can initialise the main table with the correct number of rows.
- MeasurementSetWriter correctly handles the case where the number of selected channels / frequencies is not a multiple of the averaging factors

# Version 0.3.4 (2022-01-24)

- add optional progress bars for measurement sets.

# Version 0.3.3 (2022-01-19)

- bug fixes for measurement sets when averaging

# Version 0.3.2 (2022-01-14)

- implement averaging standalone and in VisWritable
- impl Display for LatLngHeight, RADec

# Version 0.3.1 (2022-01-10)

- tweak dependency versions
- Re-export more crates
- Actually do something with cuda-static feature
- Add Cargo.lock to the gitignore
- Don't label assert lines as partially covered
- add Dockerfile

# Version 0.3.0 (2021-12-17)

- use Rust 2021 Edition
- add cuda convenience functions
- implement io feature, make rubbl optional
- import approx and ndarray directly, not rubbl's
- slightly faster MS IO

# Version 0.2.3 (2021-12-07)

- update mwalib to 0.12
- update ndarray to ">=0.15.4,<0.16"
- write to feed table
- measurement sets complete
- remove hard dependency on approx

# Version 0.2.2 (2021-11-17)

- switch from mwa_rust_core to marlu

# Version 0.2.1 (2021-11-17)

- added measurement sets
- use rubbl_casatables 0.6.0
