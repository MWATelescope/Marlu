# Marlu

<div class="bg-gray-dark" align="center" style="background-color:#24292e">
<img src="img/marlu_logo.png" alt="marlu logo" height="200px"/>
<br/>
<br/>
<a href="https://docs.rs/crate/marlu"><img src="https://docs.rs/marlu/badge.svg" alt="docs"></a>
<img src="https://github.com/MWATelescope/Marlu/workflows/Cross-platform%20tests/badge.svg" alt="Cross-platform%20tests">
<a href="https://codecov.io/gh/MWATelescope/Marlu">
  <img src="https://codecov.io/gh/MWATelescope/Marlu/branch/main/graph/badge.svg?token=CYMROMUKRI"/>
</a>
</div>

Convenience Rust code that handles coordinate transformations, Jones matrices,
etc.

## Prerequisites

- A Rust compiler with a version >= 1.60.0

  ```bash
  $ rustc -V
  rustc 1.63.0 (4b91a6ea7 2022-08-08)
  ```

  https://www.rust-lang.org/tools/install

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

To link a system-provided static library, use e.g. `CFITSIO_STATIC=1`. To link
all system-provided static libraries, use `PKG_CONFIG_ALL_STATIC=1`. To build
all C libraries and link statically (currently only `cfitsio`), use the
`all-static` feature.

## Acknowledgement

This scientific work uses data obtained from the Murchison Radio-astronomy Observatory. We
acknowledge the Wajarri Yamatji people as the traditional owners of the Observatory site.

This repo is approved by...

<img src="https://github.com/MWATelescope/Birli/raw/main/img/CIRA_Rust_Evangelism_Strike_Force.png" height="200px" alt="CIRA Rust Evangelism Strike Force logo">
