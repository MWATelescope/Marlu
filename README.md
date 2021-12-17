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
- A Rust compiler with a version >= 1.56.0

  ```bash
  $ rustc -V
  rustc 1.57.0 (f1edd0429 2021-11-29)
  ```

  https://www.rust-lang.org/tools/install

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

If using the `cuda` feature (false by default):

- [CUDA](https://docs.nvidia.com/cuda/index.html#installation-guides)
  - Ubuntu: Follow the instructions [here](https://developer.nvidia.com/cuda-downloads)
  - Arch: `cuda`
  - The library directory can be specified manually with `CUDA_LIB`
  - If not specified, `CUDA_LIBRARY_PATH` and the `/opt/cuda` and
    `/usr/local/cuda` directories are
    [searched](https://github.com/rust-cuda/cuda-sys/blob/3a973786b3482e3fdfd783cd692fbc3c665d5c11/cuda-config/src/lib.rs#L19-L46).
  - If `CUDA` is available, use `--features=cuda-static` to link it statically.

To link a system-provided static library, use e.g. `ERFA_STATIC=1`. To link all
system-provided static libraries, use `PKG_CONFIG_ALL_STATIC=1`. To build all C
libraries and link statically, use the `all-static` feature.

## Acknowledgement

This scientific work uses data obtained from the Murchison Radio-astronomy Observatory. We
acknowledge the Wajarri Yamatji people as the traditional owners of the Observatory site.

This repo is approved by...

<img src="https://github.com/MWATelescope/Birli/raw/main/img/CIRA_Rust_Evangelism_Strike_Force.png" height="200px" alt="CIRA Rust Evangelism Strike Force logo">