#!/usr/bin/env bash
# Regenerate src/pydada2/_core.pyi from the compiled extension (nanobind stubgen).
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."
uv run python -m nanobind.stubgen -m pydada2._core -o src/pydada2/_core.pyi
# nanobind emits bare `dict`/`Callable` for dict returns and optional callbacks.
sed -i '1i # mypy: disable-error-code="type-arg"' src/pydada2/_core.pyi
uv run ruff format src/pydada2/_core.pyi
