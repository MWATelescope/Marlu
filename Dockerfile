FROM ubuntu:21.04

ENV DEBIAN_FRONTEND=noninteractive
ARG DEBUG
RUN apt-get update \
    && apt-get install -y \
        aoflagger-dev \
        build-essential \
        clang \
        curl \
        git \
        jq \
        lcov \
        libcfitsio-dev \
        liberfa-dev \
        libssl-dev \
        pkg-config \
        unzip \
        zip
RUN test -z "$DEBUG" || ( \
        apt-get install -y vim gdb \
    )
RUN apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Get Rust
RUN mkdir -m755 /opt/rust /opt/cargo
ENV RUSTUP_HOME=/opt/rust CARGO_HOME=/opt/cargo PATH=/opt/cargo/bin:$PATH
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y

# Get cargo make
RUN cargo install --force cargo-make cargo-binutils

ADD . /app
WORKDIR /app

# setup the toolchain used for coverage analysis
RUN rustup toolchain install nightly-2022-01-14 --component llvm-tools-preview --profile minimal \
    && cargo +nightly-2022-01-14 install --force cargo-make --locked --version '=0.35.7' \
    && cargo +nightly-2022-01-14 install --force cargo-llvm-cov --locked --version '=0.1.15'

ENTRYPOINT [ "/bin/bash" ]