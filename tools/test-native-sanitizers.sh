#!/usr/bin/env bash
# Build and test an isolated copy with GCC ASan, LSan, and UBSan on Linux.
set -euo pipefail

dbscan_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
dbscan_work=$(mktemp -d "${TMPDIR:-/tmp}/dbscan-sanitizers.XXXXXX")
dbscan_r_home=$(R RHOME)
read -r -a dbscan_compiler <<< "$(R CMD config CXX)"
dbscan_asan=$("${dbscan_compiler[0]}" -print-file-name=libasan.so)
dbscan_cxx_runtime=$("${dbscan_compiler[0]}" -print-file-name=libstdc++.so)
if [[ ! -f "$dbscan_asan" || ! -f "$dbscan_cxx_runtime" ]]; then
  echo "This runner requires GCC's shared AddressSanitizer runtime on Linux." >&2
  exit 1
fi
echo "Sanitizer build, installation log, and library: $dbscan_work"

"${dbscan_compiler[0]}" -std=c++17 -O1 -g \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$dbscan_root/src" "$dbscan_root/tools/test-ann-points.cpp" \
  -o "$dbscan_work/test-ann-points"

mkdir -p "$dbscan_work/source" "$dbscan_work/library"
cp "$dbscan_root/DESCRIPTION" "$dbscan_root/NAMESPACE" "$dbscan_work/source/"
for dbscan_directory in R src inst man tests data; do
  cp -R "$dbscan_root/$dbscan_directory" "$dbscan_work/source/"
done

cat > "$dbscan_work/Makevars" <<'EOF'
CXXFLAGS = -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer
CXX11FLAGS = $(CXXFLAGS)
CXX14FLAGS = $(CXXFLAGS)
CXX17FLAGS = $(CXXFLAGS)
CXX20FLAGS = $(CXXFLAGS)
CXX23FLAGS = $(CXXFLAGS)
LDFLAGS = -fsanitize=address,undefined
EOF

if ! R_MAKEVARS_USER="$dbscan_work/Makevars" MAKEFLAGS="${MAKEFLAGS:--j2}" \
  R CMD INSTALL --preclean --no-test-load --no-multiarch \
  --library="$dbscan_work/library" "$dbscan_work/source" \
  > "$dbscan_work/install.log" 2>&1; then
  cat "$dbscan_work/install.log"
  exit 1
fi

# Invoke R directly: preloading ASan into the R shell launcher would also
# instrument its sed/other subprocesses and report unrelated launcher leaks.
# Keep LD_PRELOAD scoped to this process, not to installation or compilation.
# Load the C++ runtime at startup too, so ASan can intercept C++ throws.
ASAN_OPTIONS="${ASAN_OPTIONS:-detect_leaks=1:halt_on_error=1}" \
UBSAN_OPTIONS="${UBSAN_OPTIONS:-print_stacktrace=1:halt_on_error=1}" \
  "$dbscan_work/test-ann-points"

R_HOME="$dbscan_r_home" \
LD_PRELOAD="$dbscan_asan:$dbscan_cxx_runtime${LD_PRELOAD:+:$LD_PRELOAD}" \
ASAN_OPTIONS="${ASAN_OPTIONS:-detect_leaks=1:halt_on_error=1}" \
UBSAN_OPTIONS="${UBSAN_OPTIONS:-print_stacktrace=1:halt_on_error=1}" \
DBSCAN_RANDOMIZED_ITERATIONS="${DBSCAN_RANDOMIZED_ITERATIONS:-100}" \
  "$dbscan_r_home/bin/exec/R" --vanilla --slave \
  -f "$dbscan_root/tools/test-native-sanitizers.R" \
  --args "$dbscan_work/library" "$dbscan_work/source/tests/testthat"
