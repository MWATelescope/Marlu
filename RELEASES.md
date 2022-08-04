<!-- markdownlint-disable=MD025 -->

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
