#!/usr/bin/env bash
# Run the C++-facing tests against an AddressSanitizer build of pydada2._core.
# Uses a throwaway venv because the editable install's import hook would shadow the ASan build.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."
WORK="${TMPDIR:-/tmp}/pydada2-asan"
PY="$UV_PROJECT_ENVIRONMENT/bin/python"
cmake -S . -B build/asan -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_CXX_FLAGS="-fsanitize=address -fno-omit-frame-pointer" \
  -DCMAKE_SHARED_LINKER_FLAGS="-fsanitize=address" \
  -Dnanobind_DIR="$("$PY" -m nanobind --cmake_dir)" -DPython_EXECUTABLE="$PY" \
  -DCMAKE_C_COMPILER_LAUNCHER= -DCMAKE_CXX_COMPILER_LAUNCHER= >/dev/null
cmake --build build/asan >/dev/null
rm -rf "$WORK/pkg" && mkdir -p "$WORK/pkg"
cp -r src/pydada2 "$WORK/pkg/" && cp build/asan/_core*.so "$WORK/pkg/pydada2/"
if [ ! -x "$WORK/venv/bin/python" ]; then
  uv venv -q "$WORK/venv" --python "$("$PY" -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')"
fi
# Every run (a no-op when satisfied), so a cached venv picks up new dependencies.
VIRTUAL_ENV="$WORK/venv" uv pip install -q numpy pandas dnaio xopen click pytest hypothesis \
  scikit-learn scikit-misc matplotlib
PRELOAD="$(gcc -print-file-name=libasan.so) $(gcc -print-file-name=libstdc++.so)"
LOG="$WORK/lsan.log"
set +e
env LD_PRELOAD="$PRELOAD" PYTHONPATH="$WORK/pkg:$PWD/tests" ASAN_OPTIONS=detect_leaks=1 \
  "$WORK/venv/bin/python" -m pytest tests/parity tests/unit -p no:cacheprovider -o addopts="" \
  -k "not interrupt_raises and not package_version" -q "$@" >"$LOG" 2>&1
set -e
grep -E "passed|failed" "$LOG" | tail -1
if grep -q "ERROR: AddressSanitizer: [a-z-]*" "$LOG" && ! grep -q "detected memory leaks" "$LOG"; then
  echo "ASan memory errors found, see $LOG"; exit 1
fi
python3 - "$LOG" <<'PY'
import re, sys
log = open(sys.argv[1]).read()
blocks = re.split(r"\n(?=(?:Direct|Indirect) leak of)", log)
ours = [b for b in blocks if "/cpp/core/" in b or "/cpp/bindings/" in b]
print(f"leak reports involving pydada2 C++: {len(ours)}")
for b in ours[:5]:
    print(b[:2000])
sys.exit(1 if ours else 0)
PY
