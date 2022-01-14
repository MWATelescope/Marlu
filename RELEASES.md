<!-- markdownlint-disable=MD025 -->

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