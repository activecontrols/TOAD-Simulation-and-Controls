#!/usr/bin/env bash
# Verification only: host build + smoke test, flight-flags standalone build, optional Cortex-M7 compile.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POLYQP="$(cd "$HERE/.." && pwd)"
CPP="$(cd "$POLYQP/.." && pwd)"
EIGEN="$CPP/Extras/eigen-3.4.0/eigen-3.4.0"
CXX="${CXX:-c++}"
OUT="${OUT:-$HERE/build}"

if [ ! -f "$EIGEN/Eigen/Dense" ]; then
    echo "Eigen 3.4.0 not found at $EIGEN" >&2
    exit 1
fi
mkdir -p "$OUT"
cd "$POLYQP"

FLIGHT_SRCS=(PolyQP.cpp resolveWaypoints.cpp assembleQP.cpp solveLeg.cpp allocateTimes.cpp evaluateConstraints.cpp flatnessInvert.cpp
    allocateControls.cpp sampleOutputs.cpp PolyEval.cpp buildVehicleLimits.cpp TOADWaypoint.cpp ComputeJtot.cpp DCM_Quat_Conversion.cpp AbortGuidance.cpp abortFallback.cpp quatRot.cpp)

echo "== 1. Host build + smoke test (verification flags) [$CXX]"
"$CXX" -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -ffp-contract=off \
    -DEIGEN_NO_MALLOC -DPOLYQP_DEBUG_BOUNDS -DPOLYQP_ERROR_MESSAGES \
    -isystem "$EIGEN" -I. -I.. -I../Extras \
    "${FLIGHT_SRCS[@]}" verify/abort_smoke.cpp -o "$OUT/abort_smoke"
"$OUT/abort_smoke"

echo "== 2. Flight-flags build of cpp/PolyQP/*.cpp alone (no verification defines)"
for f in "${FLIGHT_SRCS[@]}"; do
    "$CXX" -std=c++17 -O2 -Wall -Wextra -Werror -ffp-contract=off -isystem "$EIGEN" -I. -I.. -I../Extras -c "$f" -o "$OUT/flight_${f%.cpp}.o"
done
echo "flight-flags build OK (${#FLIGHT_SRCS[@]} translation units)"

if command -v arm-none-eabi-g++ >/dev/null 2>&1; then
    echo "== 3. Cortex-M7 compile (no link)"
    for f in "${FLIGHT_SRCS[@]}"; do
        arm-none-eabi-g++ -std=c++17 -mcpu=cortex-m7 -mfpu=fpv5-d16 -mfloat-abi=hard -mthumb -O2 -ffp-contract=off -fno-exceptions -fno-rtti \
            -Wall -Wextra -Werror -isystem "$EIGEN" -I. -I.. -I../Extras -c "$f" -o "$OUT/m7_${f%.cpp}.o"
    done
    arm-none-eabi-size "$OUT"/m7_*.o
else
    echo "== 3. Cortex-M7 compile skipped: arm-none-eabi-g++ not found"
fi
