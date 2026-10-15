#!/usr/bin/env bash
# Runs once after the dev container is created (see .devcontainer/devcontainer.json).
set -euo pipefail

# Named volumes are created root-owned; hand them to the vscode user.
sudo chown -R vscode:vscode /home/vscode/.claude /home/vscode/.cache /home/vscode/R

# Work from the repository root regardless of where the container mounts it.
cd "$(dirname "$(readlink -f "$0")")/.."

# Git safety for bind mounts owned by a different host UID (Windows/macOS).
git config --global --add safe.directory "$PWD"

# 1. Install the original R package (reference implementation) from source, if present.
#    Before the restructure (spec.md, phase 0) it lives at the repo root; afterwards
#    in reference/r-dada2/.
if [ -f reference/r-dada2/DESCRIPTION ]; then
    R_PKG_DIR=reference/r-dada2
elif [ -f DESCRIPTION ]; then
    R_PKG_DIR=.
else
    R_PKG_DIR=""
fi
if [ -n "$R_PKG_DIR" ]; then
    echo ">> Installing reference R package from $R_PKG_DIR"
    R CMD INSTALL --no-docs --no-multiarch --library="$R_LIBS_USER" "$R_PKG_DIR" \
        || echo "!! R CMD INSTALL failed; parity fixtures cannot be regenerated until fixed."
fi

# 2. Sync the Python project (builds the nanobind extension in editable mode).
if [ -f pyproject.toml ]; then
    echo ">> uv sync"
    uv sync --all-groups || echo "!! uv sync failed; run it manually to see the error."
fi

echo ">> Done. Start Claude Code with: claude"
