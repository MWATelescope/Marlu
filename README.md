# Marlu

<!-- markdownlint-disable MD033 -->
<div class="bg-gray-dark" align="center" style="background-color:#24292e">
<img src="img/marlu_logo.png" alt="marlu logo" height="200px"/>
<br/>
<br/>
<a href="https://docs.rs/crate/marlu"><img src="https://docs.rs/marlu/badge.svg" alt="docs"></a>
<img src="https://github.com/MWATelescope/Marlu/workflows/Cross-platform%20tests/badge.svg" alt="Cross-platform%20tests">
<a href="https://codecov.io/gh/MWATelescope/Marlu">
  <img src="https://codecov.io/gh/MWATelescope/Marlu/branch/main/graph/badge.svg?token=CYMROMUKRI" alt="codecov"/>
</a>
<a href="https://crates.io/crates/marlu"><img src="https://img.shields.io/badge/rustc-1.88-orange.svg" alt="rustc"/></a>
</div>

Convenience Rust code that handles coordinate transformations, Jones matrices,
visibility averaging, and writing visibilities to uvfits files and Measurement
Sets. It is the shared base of [Birli](https://github.com/MWATelescope/Birli)
and [hyperdrive](https://github.com/MWATelescope/mwa_hyperdrive).

API documentation is on [docs.rs](https://docs.rs/marlu), and changes between
versions are listed in [RELEASES.md](RELEASES.md).

## Prerequisites

- Rust (and Cargo) version >= 1.88.0, the `rust-version` in `Cargo.toml`

```bash
cargo -V
```

<https://www.rust-lang.org/tools/install>

### Optional prerequisites

If using the `cfitsio` feature, which the `mwalib` feature (on by default)
turns on:

- [cfitsio](https://heasarc.gsfc.nasa.gov/docs/software/fitsio/)
  - Ubuntu: `libcfitsio-dev`
  - Arch: `cfitsio`
  - Library and include dirs can be specified manually with `CFITSIO_LIB` and
    `CFITSIO_INC`
  - If not specified, `pkg-config` is used to find the library.
  - Use `--features=cfitsio-static` to build the library automatically. Requires
    a C compiler and `autoconf`.

If using the `ms` feature (on by default), a C++ compiler is needed to build
the casacore tables code bundled with
[rubbl_casatables](https://crates.io/crates/rubbl_casatables).

To link a system-provided static library, use e.g. `CFITSIO_STATIC=1`. To link
all system-provided static libraries, use `PKG_CONFIG_ALL_STATIC=1`. To build
all C libraries and link statically (currently only `cfitsio`), use the
`all-static` feature.

## Features

| Feature | Default | What it provides |
|---|---|---|
| `mwalib` | yes | Conversions from [mwalib](https://github.com/MWATelescope/mwalib) contexts. Turns on `cfitsio`. |
| `ms` | yes | The Measurement Set writer. |
| `cfitsio` | through `mwalib` | The uvfits writer. Links against cfitsio. |
| `approx` | no | `approx` traits on data types. |
| `serde` | no | `Serialize` and `Deserialize` on data types. |
| `cfitsio-static` | no | Builds cfitsio from source and links it statically. |
| `all-static` | no | Builds and statically links every C library (currently only cfitsio). |

## Acknowledgement

This scientific work uses data obtained from the Murchison Radio-astronomy Observatory. We
acknowledge the Wajarri Yamatji people as the traditional owners of the Observatory site.

This repo is approved by...

<img src="https://github.com/MWATelescope/Birli/raw/main/img/CIRA_Rust_Evangelism_Strike_Force.png" height="200px" alt="CIRA Rust Evangelism Strike Force logo">
