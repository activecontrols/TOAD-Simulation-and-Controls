# PolyQP Abort Guidance: C++ (ASTRA)

This directory is a line-by-line C++ translation of the ASTRAv2 abort trajectory generator (PolyQP) from MATLAB. Open any `.cpp` next to its `.m` file: the functions, local functions, variable names, statement order and every original MATLAB comment correspond one to one.

When the vehicle aborts, the guidance builds three waypoints from the estimated state and calls `PolyQP`. That is a closed-form minimum-crackle (r = 5, degree-9) polynomial planner. It:
- works in centre-of-percussion "P-space";
- inverts differential flatness to get the 15-state trajectory and 4 feedforward inputs;
- stretches segment times until 9 actuator and kinematic limits are met.

The MATLAB is the specification. It went through three branches:
- **`abort-cpp`, the transpile itself, fixed nothing.** Anything that looked wrong was marked with a `TRANSPILE-FLAG` and listed in `TRANSPILE_FLAGS.md`.
- **`abort-fixes` fixed the problems that would break an abort in flight**, in the MATLAB first and mirrored 1:1 in the C++ (`FIX-Fxx` comments; `ABORT_FIXES.md`).
- **`abort-astra-integration` made the solver ASTRA only, NaN-free and fully static, and added the flight harness:**
  - vehicle 1 (TOAD) and every polynomial order but r = 5 were removed;
  - the constants are compiled in;
  - free boundary conditions are boolean masks instead of NaN;
  - every array is fixed-size, the grids are streamed, and all buffers share one 25 KB workspace;
  - the generic helpers moved to `cpp/matlab_helpers.cpp` / `cpp/matlab_funcs.h`;
  - the start-of-abort acceleration and jerk come from attitude and thrust;
  - the trigger was removed, because the caller decides when to abort;
  - **`cpp/Abort.h`** is what the flight software calls. See **`cpp/INTEGRATION.md`**.

## What is in this folder

### Flight code

| File | What it is |
|---|---|
| `../Abort.h`, `../Abort.cpp` | The flight harness: `Abort::plan(x_est, u_cmd, a_dist, time[], states[][15], inputs[][4], capacity)`. It converts to double, runs the abort and writes the result as float straight into the caller's trajectory arrays. It never fails: it writes either a plan or the fallback. |
| `AbortGuidance.cpp / .h` | `Guidance/AbortGuidance.m`, the abort block (formerly in `Guidance/TrajectoryManager.m`): waypoints, start boundary condition, the `PolyQP` call, status and fallback. Called by `Abort::plan`. |
| `abortFallback.cpp` | `Guidance/PolynomialQP/abortFallback.m` (FIX-F08): closed-form vertical descent, used when PolyQP fails or the caller's arrays are too small. |
| `PolyQP.cpp` | `PolyQP.m`: `PolyQP`, `solveAllLegs`, `buildReport`; it also owns the shared static workspace. |
| `resolveWaypoints.cpp` | `resolveWaypoints.m`: waypoints to knots, holds, cruise splits and legs. |
| `assembleQP.cpp` | `assembleQP.m`: normalised-time QP matrices `A0`, `Q0`, `H0` (10 × 10). |
| `solveLeg.cpp` | `solveLeg.m`: closed-form QP solve for one leg, with `R_PP` solved as an identity-padded 5 × 5 system. |
| `allocateTimes.cpp` | `allocateTimes.m`: gradient time redistribution (with `softmaxTimes`). |
| `evaluateConstraints.cpp` | `evaluateConstraints.m`: utilisation and dilation, streamed over the evaluation grid. |
| `flatnessInvert.cpp` | `flatnessInvert.m`: the per-point `flatnessInvertPoint` (used by the streamed stages), and the grid version for verification. |
| `allocateControls.cpp` | `allocateControls.m`: wrench to gimbal angles, thrust and roll. |
| `sampleOutputs.cpp` | `sampleOutputs.m`: the 0.1 s output grid, written node by node into the caller's arrays. |
| `PolyEval.cpp` | `PolyEval.m`: k-th derivative of the piecewise polynomial. |
| `buildVehicleLimits.cpp` | `buildVehicleLimits.m`: planner limits and flatness parameters (ASTRA values). |
| `TOADWaypoint.cpp` | `TOADWaypoint.m`: waypoint builder; the `arguments` defaults live in `TOADWaypointOptions`. |
| `ASTRA_PolyQP_Constants.h` | `constexpr` ASTRA constants from `LoadTOADParams.m` (m_dry, g, rTB, J, MaxThrust), with the comparison against the flight controller's values. |
| `PolyQP.h`, `PolyQP_types.h` | Declarations, structs (with the "fixed" masks), `PolyQPStatus`, `NodeSink`, the workspace. |
| `PolyQP_config.h` | Capacities (3 waypoints, 3 segments, 4 knots, 5 free unknowns, 601 output nodes) and their derivation from `ABORT_MAX_START_ALT_M`, with `static_assert`s. |
| `TRANSPILE_FLAGS.md` | Every flag and its status, the deviations, the static buffers (non-reentrancy) and the memory footprint. |
| `ABORT_FIXES.md` | The abort fixes, the decisions (D1–D10) and blocker B1, with the integration status at the top. |

These helpers moved into `cpp/` in Phase 1, in the global namespace:
- `ComputeJtot` (reduced to `J_tot = J`, `CGz = rTB`), the double `quatRot` / `DCM_Quat_Conversion` overloads, and the MATLAB-semantics helpers (`M1`, `V1`, `R1`, `Array3`, `linspace`, `cumsum`, `diff`, `findFirst`, `POLYQP_TRY`, …).
- The templates and inline functions are in `matlab_funcs.h`; the rest is in `matlab_helpers.cpp`.
- All fixed-size Eigen typedefs (`Vector15d`, `Matrix10_10d`, `Matrix15_Nd`, …) are in `matlab_funcs.h` too.

### Verification only (never part of the flight build)

| File | What it is |
|---|---|
| `verify/build_host.sh` | Host smoke test and harness test (`-Werror -DEIGEN_NO_MALLOC`), flight-flags build with per-object RAM table, optional Cortex-M7 compile, and the host solve time over the verification starts. |
| `verify/abort_smoke.cpp` | Runs `AbortGuidance` on three ASTRA states plus a capacity fallback; no MATLAB needed. |
| `../verify/abort_harness_test.cpp` | `Abort::plan` into arrays shaped like `ASTRA_Trajectory_Data.h`: all three status paths, first and last row checks, a mocked follower sequence, and 0 heap allocations. |
| `verify/abort_timing.cpp` | Times `Abort::plan` over `results/abort_cases_X_est.csv`, which `verify_polyqp` writes. |
| `verify/polyqp_mex.cpp` | MATLAB MEX dispatcher: `polyqp_mex('<FunctionName>', args...)`. It converts MATLAB NaN ↔ masks and returns `[out..., status, msg]`. |
| `verify/build_polyqp_mex.m` | Builds `polyqp_mex` (the flight capacities) or, with `'wide'`, `polyqp_mex_wide` (16 segments and a 3999 grid, for the missions and fuzzed legs). |
| `verify/verify_polyqp.m` | Per-function checks of the C++ against the `.m` files, plus every fix regression. Writes `verify/verify_polyqp_report.md`. |
| `verify/fix_regressions.m` | One regression section per fix (F07, F08, F01, F03, F27, and BC for the Phase 5 boundary condition). |
| `verify/abortWaypoints.m`, `verify/abortPair.m`, `verify/abortFallbackOnGrid.m` | The waypoint lines of `Guidance/AbortGuidance.m`, copied verbatim and asserted on every run; one abort start through MATLAB and the MEX; the fallback on a given grid. |
| `verify/check_abort_codegen.m` | Confirms `Guidance/AbortGuidance.m` still generates code (MATLAB Coder). |
| `verify/rppEqAsSolved.m` | Rebuilds the R_PP_eq that `solveLeg` factorises, for condition numbers and roundoff floors. |
| `verify/diagnose_aborts.m`, `verify/diagnosis_report.md` | The `abort-fixes` Phase-1 diagnostics, kept as a record. |
| `verify/results/`, `verify/baseline_src/` | The frozen pre-fix baseline and the pre-change sources used by the "before" comparisons (`abortWaypointsPhase4.m`, `flatnessInvert_pre.m`). |

## Building for flight

See `cpp/INTEGRATION.md` for the file list, the call and the trajectory arrays. In short:

```
-std=c++17 -O2 -ffp-contract=off -isystem <eigen> -I<cpp/PolyQP> -I<cpp> -I<cpp/Extras>
```

- **Define none** of `POLYQP_DEBUG_BOUNDS`, `POLYQP_ERROR_MESSAGES` or `EIGEN_NO_MALLOC`: they enable verification-only behaviour (index asserts, error text, allocation check). Do not override the `POLYQP_MAX_*` capacities; `Abort::MAX_POINTS` asserts against them.
- **Never** build with `-ffast-math` or `-ffinite-math-only`. The solver no longer uses NaN, but it is verified only with IEEE semantics and `-ffp-contract=off`, and no guard stops a fast-math build.
- **Memory:** 25,056 B of static RAM; the caller's 601-row abort buffer is 48 KB, and the host stack estimate is about 13 KB.
- **Not reentrant:** one static workspace. Call from one task only, never from an ISR.

## Running verification

### Host (no MATLAB)

```
cd TOAD/cpp/PolyQP/verify
./build_host.sh
```

It must print `SMOKE PASS` and `abort_harness_test: PASS`, then the RAM table, then the solve times, if `verify_polyqp` has written the case file.

### MATLAB checks

```
cd TOAD/cpp/PolyQP/verify
matlab -batch "build_polyqp_mex; build_polyqp_mex('wide'); verify_polyqp"
```

On macOS the binary is `/Applications/MATLAB_R20xxx.app/bin/matlab`. The run takes about 35 s. It compares the following with the MATLAB `.m` files:
- 200 random upright abort starts and 40 tilted, rotating starts (`PolyQP` and `AbortGuidance`, at the flight node capacity);
- every public function;
- 500 fuzzed `solveLeg` cases, about 150 fuzzed `allocateTimes` legs and three missions;
- the identity-padded `R_PP` solve against the last variable-size build (`polyqp_mex_p3`, if present). To rebuild it, run `build_polyqp_mex` in a worktree at commit `bbbf333` (Phase 3) and copy its `polyqp_mex.mex*` here as `polyqp_mex_p3.mex*`.

It then runs every `fix_regressions` section, writes `verify_polyqp_report.md`, and errors if any check fails. `check_abort_codegen` (needs MATLAB Coder) is separate.

The tolerances are the spec's, and they are not loosened. A case that misses one is accepted only if the difference is within that case's measured roundoff floor: for the linear solves, MATLAB's own error against a 60-digit `vpa` solve. The report lists every such case. The `vpa` reference requires the Symbolic Math Toolbox.

### Precision

The solver is `double` throughout. Float appears only at the `Abort::plan` boundary, as the flight types and trajectory arrays. The float solver build was removed with the helper move (`ABORT_FIXES.md` has its last measurement).

## Troubleshooting guide

| Problem | Fix |
|---|---|
| `Abort::plan` returns `FALLBACK` with `solver_status = 8` (`CapacityExceeded`) | The plan needed more than 601 output nodes (T_total > 60 s), or `capacity < Abort::MAX_POINTS`. The vehicle flies the fallback descent. If this happens inside the envelope, revisit `ABORT_MAX_START_ALT_M` and the allowance in `PolyQP_config.h`. |
| `static_assert` in `PolyQP_config.h` or `Abort.cpp` | A capacity was changed without the others. The grid must cover the derived worst case, and `Abort::MAX_POINTS` must equal `POLYQP_MAX_NODES`. |
| MATLAB and C++ differ only on the target | Check for `-ffast-math` and that `-ffp-contract=off` reaches the PolyQP sources. |
| `polyqp_mex` not found in MATLAB | Run `build_polyqp_mex` from `verify/` first; after editing C++, run `clear mex` and rebuild. |
