#!/usr/bin/env bash
# Verification only: host build + smoke test + harness test, flight-flags standalone build and RAM report, optional Cortex-M7
# compile, host solve time over the verification starts.
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
    allocateControls.cpp sampleOutputs.cpp PolyEval.cpp buildVehicleLimits.cpp TOADWaypoint.cpp AbortGuidance.cpp abortFallback.cpp ../matlab_helpers.cpp ../Abort.cpp)

echo "== 1. Host build + smoke test (verification flags) [$CXX]"
"$CXX" -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -ffp-contract=off \
    -DEIGEN_NO_MALLOC -DPOLYQP_DEBUG_BOUNDS -DPOLYQP_ERROR_MESSAGES \
    -isystem "$EIGEN" -I. -I.. -I../Extras \
    "${FLIGHT_SRCS[@]}" verify/abort_smoke.cpp -o "$OUT/abort_smoke"
"$OUT/abort_smoke"

echo "== 1b. Harness test: Abort::plan into ASTRA_Trajectory_Data.h-shaped arrays (-Werror -DEIGEN_NO_MALLOC, no verification defines)"
"$CXX" -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -ffp-contract=off -DEIGEN_NO_MALLOC \
    -isystem "$EIGEN" -I. -I.. -I../Extras \
    "${FLIGHT_SRCS[@]}" ../verify/abort_harness_test.cpp -o "$OUT/abort_harness_test"
"$OUT/abort_harness_test"

echo "== 2. Flight-flags build of cpp/PolyQP/*.cpp alone (no verification defines)"
for f in "${FLIGHT_SRCS[@]}"; do
    "$CXX" -std=c++17 -O2 -Wall -Wextra -Werror -ffp-contract=off -isystem "$EIGEN" -I. -I.. -I../Extras -c "$f" -o "$OUT/flight_$(basename "${f%.cpp}").o"
done
echo "flight-flags build OK (${#FLIGHT_SRCS[@]} translation units)"

echo "== 2b. Static RAM per object (flight-flags build): .bss + .data + common, bytes"
TOTAL=0
for f in "${FLIGHT_SRCS[@]}"; do
    o="$OUT/flight_$(basename "${f%.cpp}").o"
    if size -m "$o" >/dev/null 2>&1; then
        # Mach-O (macOS): "Section (__DATA, __bss): N"
        n=$(size -m "$o" | awk '/__DATA, __(bss|data|common)\)/ {if (match($0, /: [0-9]+/)) s += substr($0, RSTART + 2, RLENGTH - 2)} END {print s + 0}')
    else
        # ELF (Linux): Berkeley format data + bss
        n=$(size "$o" | awk 'NR == 2 {print $2 + $3}')
    fi
    printf "  %-28s %8d\n" "$(basename "$o")" "$n"
    TOTAL=$((TOTAL + n))
done
printf "  %-28s %8d  (target <= 262144)\n" "TOTAL" "$TOTAL"

if command -v arm-none-eabi-g++ >/dev/null 2>&1; then
    echo "== 3. Cortex-M7 compile (no link)"
    for f in "${FLIGHT_SRCS[@]}"; do
        arm-none-eabi-g++ -std=c++17 -mcpu=cortex-m7 -mfpu=fpv5-d16 -mfloat-abi=hard -mthumb -O2 -ffp-contract=off -fno-exceptions -fno-rtti \
            -Wall -Wextra -Werror -isystem "$EIGEN" -I. -I.. -I../Extras -c "$f" -o "$OUT/m7_$(basename "${f%.cpp}").o"
    done
    arm-none-eabi-size "$OUT"/m7_*.o
else
    echo "== 3. Cortex-M7 compile skipped: arm-none-eabi-g++ not found"
fi

CASES="$HERE/results/abort_cases_X_est.csv"
if [ -f "$CASES" ]; then
    echo "== 4. Host solve time of Abort::plan (flight flags, -O2) over the verification starts ($CASES)"
    "$CXX" -std=c++17 -O2 -Wall -Wextra -Werror -ffp-contract=off -isystem "$EIGEN" -I. -I.. -I../Extras \
        "${FLIGHT_SRCS[@]}" verify/abort_timing.cpp -o "$OUT/abort_timing"
    "$OUT/abort_timing" "$CASES"
else
    echo "== 4. Solve timing skipped: run verify_polyqp first (it writes $CASES)"
fi
