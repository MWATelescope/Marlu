ARG UBUNTU_VERSION=24.04
FROM ubuntu:${UBUNTU_VERSION}

ENV DEBIAN_FRONTEND=noninteractive
ARG DEBUG
RUN apt-get update \
    && apt-get install -y \
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
    wget \
    zip \
    && apt-get clean all && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get -y autoremove

# Get Rust
ARG RUST_VERSION=stable
ENV RUSTUP_HOME=/opt/rust CARGO_HOME=/opt/cargo PATH=/opt/cargo/bin:$PATH
RUN mkdir -m755 $RUSTUP_HOME $CARGO_HOME && \
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y \
    --profile=minimal \
    --default-toolchain=${RUST_VERSION}-$(uname -m)-unknown-linux-gnu
RUN cargo install --force cargo-llvm-cov

ADD . /app
WORKDIR /app

ENTRYPOINT [ "/bin/bash" ]