#!/usr/bin/env bash
# Standalone g2o builder for cg2o. No other shell scripts are required.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
EXTERNALS="$(realpath -m "${EXTERNALS:-$ROOT/externals}")"
BUILD_ROOT="$(realpath -m "${EXTERNAL_BUILD_ROOT:-$ROOT/build/externals}")"
INSTALL_ROOT="$(realpath -m "${EXTERNAL_INSTALL_ROOT:-$ROOT/install/externals}")"
BUILD_TYPE="${BUILD_TYPE:-Release}"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
G2O_REF="${G2O_REF:-eec325a1}"
G2O_URL="${G2O_URL:-https://github.com/RainerKuemmerle/g2o.git}"
G2O_DIR="$(realpath -m "${G2O_DIR:-$EXTERNALS/g2o}")"
G2O_BUILD_DIR="$BUILD_ROOT/g2o/$BUILD_TYPE"
G2O_INSTALL_DIR="$INSTALL_ROOT/g2o/$BUILD_TYPE"
MARCH_NATIVE="${MARCH_NATIVE:-ON}"
G2O_OPENMP="${G2O_OPENMP:-ON}"

fail() { echo "ERROR: $*" >&2; exit 1; }
config_file() { find "$1" -type f -name g2oConfig.cmake -print -quit 2>/dev/null || true; }
solvers_installed() {
  local solver
  for solver in solver_cholmod solver_csparse; do
    find "$G2O_INSTALL_DIR" -type f -name "libg2o_${solver}*" -print -quit 2>/dev/null | grep -q . || return 1
  done
}

CLEAN=0; REBUILD=0
for arg in "$@"; do
  case "$arg" in
    --clean) CLEAN=1; REBUILD=1 ;;
    --rebuild) REBUILD=1 ;;
    -h|--help)
      echo 'Usage: BUILD_TYPE=Debug|Release bash scripts/build_g2o.sh [--rebuild|--clean]'
      exit 0 ;;
    *) fail "Unknown argument: $arg" ;;
  esac
done
case "$BUILD_TYPE" in Debug|Release|RelWithDebInfo|MinSizeRel) ;; *) fail "Invalid BUILD_TYPE: $BUILD_TYPE" ;; esac
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || fail 'JOBS must be a positive integer'
[[ "$MARCH_NATIVE" == ON || "$MARCH_NATIVE" == OFF ]] || fail 'MARCH_NATIVE must be ON or OFF'
case "$INSTALL_ROOT" in /usr|/usr/*) fail "Refusing system install prefix: $INSTALL_ROOT" ;; esac
for cmd in cmake git sha256sum; do command -v "$cmd" >/dev/null 2>&1 || fail "Missing required command: $cmd"; done

GEN_ARGS=()
if [[ -n "${CMAKE_GENERATOR:-}" ]]; then
  GEN_ARGS=(-G "$CMAKE_GENERATOR")
elif command -v ninja >/dev/null 2>&1; then
  GEN_ARGS=(-G Ninja)
fi
ARCH_FLAGS=''
[[ "$MARCH_NATIVE" != ON ]] || ARCH_FLAGS='-march=native'
CMAKE_ARGS=(
  "${GEN_ARGS[@]}"
  -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
  -DCMAKE_INSTALL_PREFIX="$G2O_INSTALL_DIR"
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON
  "-DCMAKE_CXX_FLAGS=$ARCH_FLAGS"
  -DBUILD_WITH_MARCH_NATIVE="$MARCH_NATIVE"
  -DG2O_NO_IMPLICIT_OWNERSHIP_OF_OBJECTS=ON
  -DG2O_USE_OPENMP="$G2O_OPENMP"
  -DG2O_BUILD_APPS=OFF -DG2O_BUILD_EXAMPLES=OFF -DG2O_USE_OPENGL=OFF
  -DG2O_USE_CHOLMOD=ON -DG2O_USE_CSPARSE=ON
  -DG2O_FAST_MATH=OFF -DBUILD_UNITTESTS=OFF
)
HASH="$(printf '%s\n' "${CC:-}" "${CXX:-}" "source=$G2O_DIR" "ref=$G2O_REF" "${CMAKE_ARGS[@]}" | sha256sum | cut -c1-16)"
MARKER="$G2O_INSTALL_DIR/.macc_install_complete"

echo "=== g2o $BUILD_TYPE ==="
if [[ "$CLEAN" == 0 && "$REBUILD" == 0 ]] \
    && [[ -n "$(config_file "$G2O_INSTALL_DIR")" ]] \
    && solvers_installed \
    && [[ -f "$MARKER" ]] \
    && grep -Fqx "source=$G2O_DIR" "$MARKER" \
    && grep -Fqx "config=$HASH" "$MARKER"; then
  echo "Already installed: $G2O_INSTALL_DIR"
  exit 0
fi

# Preserve existing checkouts; pin only newly cloned sources.
if [[ -f "$G2O_DIR/CMakeLists.txt" ]]; then
  echo "Using existing g2o source: $G2O_DIR"
elif [[ -e "$G2O_DIR" ]]; then
  fail "Source path exists but is not a g2o checkout: $G2O_DIR"
else
  mkdir -p "$(dirname "$G2O_DIR")"
  git clone --recursive "$G2O_URL" "$G2O_DIR"
  [[ -z "$G2O_REF" ]] || git -C "$G2O_DIR" checkout "$G2O_REF"
  git -C "$G2O_DIR" submodule update --init --recursive
fi

# Legacy cg2o compatibility adjustment, as in the provided scripts.
HEADER="$G2O_DIR/g2o/core/sparse_optimizer.h"
if [[ -f "$HEADER" ]] && grep -Fqx '  void update(const double* update);' "$HEADER"; then
  sed -i 's/^  void update(const double\* update);/  virtual void update(const double* update);/' "$HEADER"
fi

if [[ "$CLEAN" == 1 ]]; then
  rm -rf -- "$G2O_BUILD_DIR" "$G2O_INSTALL_DIR"
fi
PREVIOUS="$(cat "$G2O_BUILD_DIR/.macc_config" 2>/dev/null || true)"
if [[ -d "$G2O_BUILD_DIR" && "$PREVIOUS" != "$HASH" ]]; then
  rm -rf -- "$G2O_BUILD_DIR"
fi
mkdir -p "$G2O_BUILD_DIR"
printf '%s\n' "$HASH" > "$G2O_BUILD_DIR/.macc_config"

cmake -S "$G2O_DIR" -B "$G2O_BUILD_DIR" "${CMAKE_ARGS[@]}"
cmake --build "$G2O_BUILD_DIR" --parallel "$JOBS"
rm -rf -- "$G2O_INSTALL_DIR"
cmake --install "$G2O_BUILD_DIR"
[[ -n "$(config_file "$G2O_INSTALL_DIR")" ]] || fail "g2oConfig.cmake not installed under $G2O_INSTALL_DIR"
solvers_installed || fail 'g2o CHOLMOD/CSparse libraries missing; verify SuiteSparse and the g2o revision'
printf 'source=%s\nconfig=%s\nbuild_type=%s\nmarch_native=%s\n' \
  "$G2O_DIR" "$HASH" "$BUILD_TYPE" "$MARCH_NATIVE" > "$MARKER"
echo "g2o ready: $G2O_INSTALL_DIR"

