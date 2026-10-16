# Development image for pyDADA2 (C++ + Python, uv + nanobind).
#
# Contents:
#   - Ubuntu 24.04 (noble) base
#   - C/C++ toolchain: gcc/g++, clang, clangd, cmake, ninja, ccache, gdb, valgrind, gcovr/lcov
#   - Python 3.13 managed by uv (with headers, needed for nanobind)
#   - optional (INSTALL_R=true): R (current CRAN release) + R dada2's dependencies as binaries
#     via r2u. Only needed to regenerate the parity fixtures or to benchmark against R
#     (scripts/install_r_dada2.sh then installs dada2 itself).
#   - Claude Code (native installer)
#
# Ubuntu is used instead of python:*-slim because r2u (binary CRAN/Bioconductor
# packages via apt) only targets Ubuntu. Python itself comes from uv.

FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
ARG PYTHON_VERSION=3.13
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=1000
ARG INSTALL_R=false

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    TZ=Etc/UTC

# --- Base system + C/C++ toolchain ------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl wget git gnupg sudo less nano jq unzip xz-utils \
        ripgrep fd-find procps openssh-client \
        build-essential gcc g++ gfortran make cmake ninja-build pkg-config ccache \
        clang clangd clang-format clang-tidy lldb \
        gdb valgrind lcov gcovr \
        zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev \
        libboost-dev \
    && rm -rf /var/lib/apt/lists/*

# --- Optional: R via CRAN apt + r2u (binary CRAN & Bioconductor packages) ----------
# Follows https://github.com/eddelbuettel/r2u/blob/master/inst/scripts/add_cranapt_noble.sh
RUN if [ "$INSTALL_R" != "true" ]; then exit 0; fi \
    && gpg --homedir /tmp --no-default-keyring --keyring /usr/share/keyrings/r2u.gpg \
        --keyserver keyserver.ubuntu.com \
        --recv-keys A1489FE2AB99A21A 67C2D66C4B1D4339 51716619E084DAB9 \
    && printf 'Types: deb\nURIs: https://r2u.stat.illinois.edu/ubuntu\nSuites: noble\nComponents: main\nArch: amd64, arm64\nSigned-By: /usr/share/keyrings/r2u.gpg\n' \
        > /etc/apt/sources.list.d/r2u.sources \
    && printf 'Types: deb\nURIs: https://cloud.r-project.org/bin/linux/ubuntu\nSuites: noble-cran40/\nComponents:\nArch: amd64, arm64\nSigned-By: /usr/share/keyrings/r2u.gpg\n' \
        > /etc/apt/sources.list.d/cran.sources \
    && printf 'Package: *\nPin: release o=CRAN-Apt Project\nPin: release l=CRAN-Apt Packages\nPin-Priority: 700\n' \
        > /etc/apt/preferences.d/99cranapt \
    && apt-get update && apt-get install -y --no-install-recommends \
        r-base-core r-base-dev \
        r-cran-rcpp r-cran-rcppparallel r-cran-ggplot2 r-cran-reshape2 \
        r-cran-jsonlite r-cran-knitr r-cran-rmarkdown \
        r-bioc-biocgenerics r-bioc-biostrings r-bioc-iranges r-bioc-xvector r-bioc-shortread \
    && rm -rf /var/lib/apt/lists/* /tmp/*

# --- uv (pinned) ------------------------------------------------------------------
COPY --from=ghcr.io/astral-sh/uv:0.12.8 /uv /uvx /usr/local/bin/

# Shared, uv-managed Python (includes headers for building C++ extensions)
ENV UV_PYTHON_INSTALL_DIR=/opt/uv/python \
    UV_PYTHON_PREFERENCE=only-managed \
    UV_PYTHON=${PYTHON_VERSION}
RUN uv python install ${PYTHON_VERSION} \
    && chmod -R a+rX /opt/uv

# --- Non-root user ----------------------------------------------------------------
# ubuntu:24.04 ships a default 'ubuntu' user with UID 1000; remove it first.
RUN if id -u ubuntu >/dev/null 2>&1; then userdel -r ubuntu; fi \
    && groupadd --gid ${USER_GID} ${USERNAME} \
    && useradd --uid ${USER_UID} --gid ${USER_GID} -m -s /bin/bash ${USERNAME} \
    && echo "${USERNAME} ALL=(root) NOPASSWD:ALL" > /etc/sudoers.d/${USERNAME} \
    && chmod 0440 /etc/sudoers.d/${USERNAME} \
    && mkdir -p /workspace && chown ${USERNAME}:${USERNAME} /workspace

USER ${USERNAME}

# --- Claude Code (native installer, as the vscode user) ----------------------------
RUN curl -fsSL https://claude.ai/install.sh | bash

# Per-user dirs (also volume mount points; pre-created so they get correct ownership)
RUN mkdir -p /home/${USERNAME}/.claude /home/${USERNAME}/.cache/uv /home/${USERNAME}/.cache/ccache \
             /home/${USERNAME}/R/library

ENV PATH="/home/${USERNAME}/.venv/bin:/home/${USERNAME}/.local/bin:${PATH}" \
    UV_LINK_MODE=copy \
    UV_PROJECT_ENVIRONMENT=/home/${USERNAME}/.venv \
    UV_CACHE_DIR=/home/${USERNAME}/.cache/uv \
    CCACHE_DIR=/home/${USERNAME}/.cache/ccache \
    CMAKE_GENERATOR=Ninja \
    CMAKE_C_COMPILER_LAUNCHER=ccache \
    CMAKE_CXX_COMPILER_LAUNCHER=ccache \
    R_LIBS_USER=/home/${USERNAME}/R/library

WORKDIR /workspace
