# PolyQP Abort Guidance: C++ Transpile

This directory is a line-by-line C++ translation of the ASTRAv2 abort trajectory generator (PolyQP) from MATLAB, at commit `5ac6bd23`. Open any `.cpp` next to its `.m` file: the functions, local functions, variable names, statement order and every original MATLAB comment correspond one to one.

When the abort fires, the guidance builds three waypoints from the estimated state and calls `PolyQP`. That is a closed-form minimum-crackle (r = 5, degree-9) polynomial planner that works in centre-of-percussion "P-space", inverts differential flatness to get the 15-state trajectory and 4 feedforward inputs, and stretches segment times until 9 actuator and kinematic limits are met.

The MATLAB is the specification. The transpile itself (branch `abort-cpp`) fixed nothing; anything that looked wrong was marked with a `TRANSPILE-FLAG` and listed in `TRANSPILE_FLAGS.md`. Branch `abort-fixes` then fixed the problems that would break an abort in flight, in the MATLAB first and mirrored 1:1 in the C++. Every changed line carries a matching `FIX-Fxx` comment; see `ABORT_FIXES.md` for each fix, the decisions (D1–D8) and the controller-side blocker B1.

## What is in this folder

### Flight code

| File | What it is |
|---|---|
| `AbortGuidance.cpp / .h` | The abort block of `Guidance/TrajectoryManager.m` (waypoints, `PolyQP` call, fallback, `Time`/`States`/`Inputs`). The persistents are explicit state in `AbortGuidanceState`; the outcome is `AbortStatus`. |
| `abortFallback.cpp` | `Guidance/PolynomialQP/abortFallback.m` (FIX-F08): closed-form vertical descent used when PolyQP fails. |
| `quatRot.cpp` | `Helper/quatRot.m` in `double` (FIX-F03: b3 = C_IB e3 from the estimated quaternion). |
| `PolyQP.cpp` | `PolyQP.m`: `PolyQP`, `solveAllLegs`, `buildReport`. |
| `resolveWaypoints.cpp` | `resolveWaypoints.m`: waypoints to knots, holds, cruise splits and legs. |
| `assembleQP.cpp` | `assembleQP.m`: normalised-time QP matrices `A0`, `Q0`, `H0`. |
| `solveLeg.cpp` | `solveLeg.m`: closed-form QP solve for one leg. |
| `allocateTimes.cpp` | `allocateTimes.m`: gradient time redistribution (with `softmaxTimes`). |
| `evaluateConstraints.cpp` | `evaluateConstraints.m`: dense-grid utilisation and dilation (with the two thrust-stretch roots). |
| `flatnessInvert.cpp` | `flatnessInvert.m`: flat outputs to states and wrench. |
| `allocateControls.cpp` | `allocateControls.m`: wrench to gimbal angles, thrust and roll. |
| `sampleOutputs.cpp` | `sampleOutputs.m`: the output grid at `DtSample`. |
| `PolyEval.cpp` | `PolyEval.m`: k-th derivative of the piecewise polynomial. |
| `buildVehicleLimits.cpp` | `buildVehicleLimits.m`: planner limits and flatness parameters. |
| `TOADWaypoint.cpp` | `TOADWaypoint.m`: waypoint builder; the `arguments` defaults live in `TOADWaypointOptions`. |
| `ComputeJtot.cpp` | `Helper/ComputeJtot.m`: inertia and CG versus propellant mass. |
| `DCM_Quat_Conversion.cpp` | `Helper/DCM_Quat_Conversion.m` in `double`, inside `namespace PolyQP`. The `float` version in `cpp/matlab_helpers.cpp` is unchanged (F15). |
| `PolyQP.h` | Declarations of the public functions. |
| `PolyQP_types.h` | `Real` (`double` by default), fixed-capacity Eigen types, `Array3`, all structs, and `PolyQPStatus`. |
| `PolyQP_config.h` | Capacities (`POLYQP_MAX_WP`, `_SEG`, `_KNOTS`, `_R`, `_NC`, `_GRID`). They default to MATLAB's values and can each be overridden with `-D`. |
| `matlab_compat.h` | Header-only MATLAB-semantics helpers: 1-based access (`M1`, `V1`, `R1`), NaN-skipping `mmax`/`mmin`, MATLAB `linspace`, `find`, `cumsum`, `diff`, `any`/`all`, plus the fast-math guard. |
| `TRANSPILE_FLAGS.md` | Every flag and its status after the fixes, the deviations, the static buffers (non-reentrancy) and the memory footprint. |
| `ABORT_FIXES.md` | The abort fixes: cause, change (file:line), acceptance result, decisions, blocker, before/after table. |

### Verification only (`verify/`, never part of the flight build)

| File | What it is |
|---|---|
| `verify/build_host.sh` | Host build plus smoke test, flight-flags standalone build, and an optional Cortex-M7 compile. |
| `verify/abort_smoke.cpp` | Runs `AbortGuidance` on 3 ASTRA states; no MATLAB needed. |
| `verify/polyqp_mex.cpp` | MATLAB MEX dispatcher: `polyqp_mex('<FunctionName>', args...)` for every public function. Returns `[out..., status, msg]`. |
| `verify/build_polyqp_mex.m` | Builds `polyqp_mex` (double) or, with `'float'`, `polyqp_mex_float` (`-DPOLYQP_REAL=float`). |
| `verify/verify_polyqp.m` | Per-function checks of the C++ against the `.m` files. Writes `verify/verify_polyqp_report.md`. |
| `verify/verify_polyqp_report.md` | The latest verification report (transpile checks plus every fix regression). |
| `verify/fix_regressions.m` | One regression section per FIX (acceptance criteria of `ABORT_FIXES.md`); run by `verify_polyqp`. |
| `verify/diagnose_aborts.m`, `verify/diagnosis_report.md` | Phase-1 abort diagnostics (infeasibility, conditioning, roll units, t = 0 steps, P-offset), before and after the fixes. |
| `verify/check_abort_codegen.m` | Confirms `TrajectoryManager.m` (with the abort branch enabled) still generates code, as the Simulink block requires. |
| `verify/matlabAbortBlock.m`, `verify/abortWaypoints.m` | Verbatim copies of the abort block and its waypoint lines (asserted against `TrajectoryManager.m` on every run). |
| `verify/rppEqAsSolved.m` | Rebuilds the R_PP_eq that `solveLeg` factorises, for condition numbers and roundoff floors. |
| `verify/results/` | The frozen pre-fix baseline (`baseline_report.md`, `baseline_cases_slim.mat`) and the diagnostics data. |
| `verify/baseline_src/flatnessInvert_pre.m` | The pre-fix `flatnessInvert` (renamed), for the FIX-F27 "unchanged" check. |

## Building for flight

Compile the 17 `.cpp` files in this folder alongside the existing flight code. They need only Eigen 3.4 and `cpp/matlab_funcs.h` (for the `Vector15` boundary type):

```
-std=c++17 -O2 -ffp-contract=off -isystem <eigen> -I<cpp/PolyQP> -I<cpp> -I<cpp/Extras>
```

- **Define none** of `POLYQP_DEBUG_BOUNDS`, `POLYQP_ERROR_MESSAGES`, `EIGEN_NO_MALLOC` or `POLYQP_REAL`. They enable verification-only behaviour (index asserts, error text, allocation check, precision override).
- **Never** build with `-ffast-math` or `-ffinite-math-only`. PolyQP uses NaN as the "free" boundary-condition marker and `inf` limits; `matlab_compat.h` stops the build if fast math is on.
- Include PolyQP headers before any other Eigen header in a translation unit (F25).
- **Memory: the default capacities do not fit a Teensy 4.x / STM32H7.** At MATLAB's capacities PolyQP needs 13.0 MB of static RAM (9.6 MB before FIX-F07); `-DPOLYQP_MAX_SEG=4 -DPOLYQP_MAX_GRID=1200` still needs 1.8 MB. Depleted-tank TOAD aborts can exceed 400 s (FIX-F29), which matters for the grid capacity. See `TRANSPILE_FLAGS.md` → Memory footprint. Choosing a flight profile (and checking it against the longest abort) is a flight-integration task.
- PolyQP is **not reentrant** (large buffers are `static`). Call it from one task only, never from an ISR.

### Calling the abort guidance

```cpp
#include "AbortGuidance.h"   // before other Eigen includes

PolyQP::AbortGuidanceState abortState;          // abortCalculated = false, tAbort = 0
PolyQP::PolyQPConstants    C = /* fields of LoadTOADParams: Vehicle, m_dry, g, rTB, J, MaxThrust, ... */;
PolyQP::AbortTimeF   Time;                      // N x 1
PolyQP::AbortStatesF States;                    // N x 15
PolyQP::AbortInputsF Inputs;                    // N x 4 (last row duplicated, as in MATLAB)
PolyQP::ReportT      Report;

PolyQP::AbortStatus st = PolyQP::AbortGuidance(t, X_est, C, abortState, Time, States, Inputs, Report /*, ABORT = 0, ABORT_PIN_THRUST = false */);
```

`ABORT` defaults to 0, the MATLAB value, so the abort never fires unless the caller passes a trigger time (F09, F31). When it fires, `Time`/`States`/`Inputs` are written **once** and `abortState.abortCalculated`/`tAbort` are set; later calls change nothing (FIX-F32). The return value (also `abortState.abortStatus`) says what was written:

- `OK`: the PolyQP plan, feasible.
- `INFEASIBLE_FLOWN`: the PolyQP plan with `Report.Feasible = 0`, flown anyway (D2); log `Report`.
- `FALLBACK`: PolyQP failed (`abortState.planStatus` holds the error, e.g. `SingularRPP` or `CapacityExceeded`) or the estimate was non-finite. The trajectory is the closed-form vertical descent of `abortFallback` (D1).

The plan uses the estimated attitude for the knot-1 P-space BC (FIX-F03) and the estimated propellant `X_est(14:15)` (FIX-F29). `ABORT_PIN_THRUST` (default false, D4) also pins the knot-1 acceleration and jerk to hover along the current thrust axis. A `Real` (`double`) overload with the same arguments exists for verification.

## Running verification

### Host smoke test (no MATLAB)

```
cd TOAD/cpp/PolyQP/verify
./build_host.sh
```

This builds the flight sources plus `abort_smoke.cpp` with `-Wall -Wextra -Wpedantic -Werror -ffp-contract=off -DEIGEN_NO_MALLOC -DPOLYQP_DEBUG_BOUNDS -DPOLYQP_ERROR_MESSAGES` and runs it (it must print `SMOKE PASS`). It then compiles each flight file alone with the flight flags, and compiles for Cortex-M7 if `arm-none-eabi-g++` is installed.

### MATLAB checks

```
cd TOAD/cpp/PolyQP/verify
matlab -batch "build_polyqp_mex; build_polyqp_mex('float'); verify_polyqp"
```

On macOS the binary is `/Applications/MATLAB_R20xxx.app/bin/matlab`. The transpile checks take about a minute and the fix regressions about ten more. It compares the following with the MATLAB `.m` files:
- 200 random abort states per vehicle plus 40 tilted, rotating, depleted-tank states per vehicle (`PolyQP` and `AbortGuidance`);
- every public function;
- 500 fuzzed `solveLeg` cases, 150 fuzzed `allocateTimes` legs and three missions.

It then runs every `fix_regressions` section and measures the float build for comparison. It writes `verify_polyqp_report.md` and errors if any double-precision check or fix regression fails. `check_abort_codegen` (needs MATLAB Coder) and `diagnose_aborts('fixed')` are separate.

The tolerances are the spec's, and they are not loosened. A case that misses one is accepted only if the difference is within that case's measured roundoff floor: for the linear solves, MATLAB's own error against a 60-digit `vpa` solve, or the textbook cond·eps bound. The report lists every such case. Requires the Symbolic Math Toolbox for the `vpa` reference.

### Precision

`Real` is `double`. Before FIX-F07 the float build (`-DPOLYQP_REAL=float`) failed with `SingularRPP` in 10 of 200 ASTRA abort cases and showed up to 26 m position error against MATLAB. After FIX-F07 it has no status failures and at most 0.021 m position error (`verify_polyqp_report.md`). The solver internals still stay in double, because the measured float error is informational only. Only the `AbortGuidance` boundary uses the flight `float` types.

## Troubleshooting guide

| Problem | Fix |
|---|---|
| `PolyQP requires IEEE NaN/Inf semantics` | Remove `-ffast-math` / `-ffinite-math-only` from the build. |
| `Include PolyQP headers before any Eigen header` | Move `#include "PolyQP.h"` / `"AbortGuidance.h"` above other includes that pull in Eigen. |
| Static assertion about object size from Eigen | Same cause: Eigen was configured before `PolyQP_config.h` set `EIGEN_STACK_ALLOCATION_LIMIT`. |
| Linker reports RAM overflow on the target | Expected at the default capacities; set a smaller profile with `-DPOLYQP_MAX_SEG=... -DPOLYQP_MAX_GRID=...` (see Memory above). |
| `AbortGuidance` returns `FALLBACK` with `planStatus = CapacityExceeded` | The trajectory needed more segments or grid points than the compiled capacities; the vehicle flies the fallback descent. Raise the capacity profile. |
| `polyqp_mex` not found in MATLAB | Run `build_polyqp_mex` from `verify/` first; after editing C++, run `clear mex` and rebuild. |
