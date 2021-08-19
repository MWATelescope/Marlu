# mwa_rust_core

<div class="bg-gray-dark" align="center" style="background-color:#24292e">
<br/>
<a href="https://docs.rs/crate/mwa_rust_core"><img src="https://docs.rs/mwa_rust_core/badge.svg" alt="docs"></a>
<img src="https://github.com/MWATelescope/giant-squid/workflows/Cross-platform%20tests/badge.svg" alt="Cross-platform%20tests">
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
