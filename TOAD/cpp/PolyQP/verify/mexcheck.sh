#!/usr/bin/env bash
# Verification only: fast syntax check of polyqp_mex.cpp against MATLAB's extern headers (no MATLAB session needed).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CPP="$(cd "$HERE/../.." && pwd)"
ML="${MATLAB_ROOT:-$(ls -d /Applications/MATLAB_R20*.app 2>/dev/null | tail -1)}"
c++ -std=c++17 -fsyntax-only -Wall -Wextra -ffp-contract=off -DPOLYQP_DEBUG_BOUNDS -DPOLYQP_ERROR_MESSAGES -D_USE_MATH_DEFINES -DMATLAB_DEFAULT_RELEASE=R2017b \
    -isystem "$ML/extern/include" -isystem "$CPP/Extras/eigen-3.4.0/eigen-3.4.0" -I"$CPP/PolyQP" -I"$CPP" -I"$CPP/Extras" "$@" "$HERE/polyqp_mex.cpp"
echo "polyqp_mex.cpp syntax OK $*"
