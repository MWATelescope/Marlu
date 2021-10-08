# mwa_rust_core

<div class="bg-gray-dark" align="center" style="background-color:#24292e">
<br/>
<a href="https://docs.rs/crate/mwa_rust_core"><img src="https://docs.rs/mwa_rust_core/badge.svg" alt="docs"></a>
<img src="https://github.com/MWATelescope/mwa_rust_core/workflows/Cross-platform%20tests/badge.svg" alt="Cross-platform%20tests">
</div>

Convenience Rust code that handles coordinate transformations, Jones matrices,
etc.

## Prerequisites
- A Rust compiler with a version >= 1.50.0

  `https://www.rust-lang.org/tools/install`

- [ERFA](https://github.com/liberfa/erfa)
  - Ubuntu: `liberfa-dev`
  - Arch: AUR package `erfa`
  - The library dir can be specified manually with `ERFA_LIB`
  - If not specified, `pkg-config` is used to find the library.
  - Use `--features=erfa-static` to build the library automatically. Requires a
    C compiler and `autoconf`.

### Optional prerequisites
If using the `mwalib` feature (true by default):

- [cfitsio](https://heasarc.gsfc.nasa.gov/docs/software/fitsio/)
  - Ubuntu: `libcfitsio-dev`
  - Arch: `cfitsio`
  - Library and include dirs can be specified manually with `CFITSIO_LIB` and
    `CFITSIO_INC`
  - If not specified, `pkg-config` is used to find the library.
  - Use `--features=cfitsio-static` to build the library automatically. Requires
    a C compiler and `autoconf`.

To link a system-provided static library, use e.g. `ERFA_STATIC=1`. To link all
system-provided static libraries, use `PKG_CONFIG_ALL_STATIC=1`. To build all C
libraries and link statically, use the `all-static` feature.

## Troubleshooting

### the trait bound `Jones<f32>: AbsDiffEq<_>` is not satisfied

if you see an error that looks like this:

```txt
error[E0277]: the trait bound `Jones<f32>: AbsDiffEq<_>` is not satisfied
     |
1029 | /         assert_abs_diff_eq!(
1030 | |             *jones_array.get((3, 3, 1)).unwrap(),
1031 | |             &Jones::from([
1032 | |                 Complex::new(rot_1_xx_3_3_re, rot_1_xx_3_3_im),
...    |
1036 | |             ])
1037 | |         );
     | |__________^ the trait `AbsDiffEq<_>` is not implemented for `Jones<f32>`
     |
     = note: this error originates in the macro `abs_diff_eq` (in Nightly builds, run with -Z macro-backtrace for more info)
```

try

```bash
cargo update
cargo update -p approx:0.5.0 --precise 0.4.0
cargo update -p ndarray:0.15.3 --precise 0.14.0
```

and add `--locked` to any cargo commands that might perform a `cargo update` (e.g. `cargo install`)