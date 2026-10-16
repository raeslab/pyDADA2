#!/usr/bin/env bash
# Runs once after the dev container is created (see .devcontainer/devcontainer.json).
set -euo pipefail

# Named volumes are created root-owned; hand them to the vscode user.
sudo chown -R vscode:vscode /home/vscode/.claude /home/vscode/.cache /home/vscode/R

# Work from the repository root regardless of where the container mounts it.
cd "$(dirname "$(readlink -f "$0")")/.."

# Git safety for bind mounts owned by a different host UID (Windows/macOS).
git config --global --add safe.directory "$PWD"

# 1. Optional: R dada2 1.36 for regenerating fixtures / R benchmarks (image built with INSTALL_R=true).
if command -v R >/dev/null 2>&1 && ! Rscript -e 'library(dada2)' >/dev/null 2>&1; then
    echo ">> Installing R dada2 1.36 (reference for the parity fixtures)"
    bash scripts/install_r_dada2.sh \
        || echo "!! R dada2 install failed; parity fixtures cannot be regenerated until fixed."
fi

# 2. Sync the Python project (builds the nanobind extension in editable mode).
echo ">> uv sync"
uv sync --all-groups || echo "!! uv sync failed; run it manually to see the error."

echo ">> Done."
