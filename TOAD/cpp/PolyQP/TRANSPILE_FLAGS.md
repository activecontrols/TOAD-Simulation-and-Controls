# PolyQP transpile: flag log

Every `TRANSPILE-FLAG` in `cpp/PolyQP/` is listed here. A flag marks something in the MATLAB that looks wrong, or a place where the C++ had to differ from the MATLAB beyond the permitted deviations (spec Section 5.8).

On branch `abort-cpp` nothing flagged was fixed: the C++ did what the MATLAB at `5ac6bd23` does. On branch `abort-fixes` the safety-relevant flags were fixed **in the MATLAB first and mirrored in the C++** (see `ABORT_FIXES.md`). A fixed flag's comment reads `TRANSPILE-FLAG Fxx (fixed by FIX-Fxx)`. Every changed line carries a `FIX-Fxx` comment in both languages. The **Status** column below gives the state on `abort-fixes`.

Branch `abort-astra-integration` then made the solver ASTRA only, NaN-free and fully fixed-size, moved the generic helpers into `cpp/matlab_helpers.cpp` / `cpp/matlab_funcs.h`, replaced F02's option with an always-on start boundary condition, removed the trigger, and added the flight harness `cpp/Abort.h` (see `README.md` and `cpp/INTEGRATION.md`). The **Status** column includes those changes; flags whose code is gone are marked *removed* or *moot*.

MATLAB paths are relative to `TOAD/`; C++ paths to `TOAD/cpp/PolyQP/`. Line numbers are at `abort-fixes` and have moved since; find a flag in the C++ with `grep -rn "TRANSPILE-FLAG Fxx"`. The abort block, formerly in `Guidance/TrajectoryManager.m`, is now `Guidance/AbortGuidance.m` (commit fd08005 replaced `TrajectoryManager.m`). IDs F17–F23 are unused: new flags were numbered from F24.

Severity: **info** (no effect on outputs), **behaviour** (changes what the planner produces compared with what was presumably intended), **safety** (can affect the vehicle in an abort).

## Pre-identified flags (F01–F16)

| ID | File:line (.m) | File:line (.cpp) | Description | Severity | Status |
| --- | --- | --- | --- | --- | --- |
| F01 | `Guidance/TrajectoryManager.m:15`, `:41` | `AbortGuidance.cpp:14`, `:48` | Abort waypoint 2 altitude is 2 m; the lead's spec says 1 m above the pad. | behaviour | **Open, DECISION-D5**: now the named constant `ABORT_HOVER_ALT = 2` (no value change). |
| F02 | `Guidance/TrajectoryManager.m:16`, `:37` | `AbortGuidance.cpp:42` | By default only position and velocity are pinned at knot 1. Accel/jerk/snap are free, so thrust and attitude feedforward can step at the abort instant (measured: up to 25° attitude, 41 % thrust). | safety | **Fixed by FIX-F02 (ASTRA integration Phase 5)**: knot 1 always pins a0 = C_IB T_B / m + g + a_dist and j0 = (w_I × C_IB T_B) / m from attitude, thrust command and disturbance (D6: Ṫ = 0); `ABORT_PIN_THRUST` is removed. As D4 predicted, most tilted starts are now `INFEASIBLE_FLOWN` (13 / 100 feasible). |
| F03 | `Guidance/PolynomialQP/PolyQP.m:140`, `:223` | `PolyQP.cpp:169`, `:270` | Knot-1 P-offset assumed b3 = e3 and was never refined. The velocity BC was applied to dP/dt, not v_cg + rho·db3/dt. | behaviour | **Fixed by FIX-F03** (`B3BC`/`B3DotBC`, abort passes the estimated attitude). |
| F04 | `Guidance/PolynomialQP/PolyQP.m:186` | `PolyQP.cpp:228` | `nPasses = 3`: for the abort, passes 2–3 recompute pass 1 exactly. Wasted work in flight. | info | Open (out of scope). |
| F05 | `Guidance/PolynomialQP/PolyQP.m:37` | `PolyQP.cpp:62` | `legCanContract` is false for every abort, so the abort never contracts below its initial time guess. | behaviour | Open (out of scope). |
| F06 | `Guidance/PolynomialQP/resolveWaypoints.m:90`, `:96` | `resolveWaypoints.cpp:108`, `:116` | `all(~isnan(dp))` is false with NaN lateral positions, so `StraightLine` on abort waypoint 3 is inert and the abort never takes the cruise split. | behaviour | **Consequence addressed by FIX-F07**: the free lateral landing position no longer degrades conditioning. `StraightLine` itself is still inert. |
| F07 | `Guidance/PolynomialQP/solveLeg.m:94` | `solveLeg.cpp:138` | Lateral R_PP was near-singular for unbalanced segment times (cond ≈ (T1/T2)^9). `SingularRPP` tripped in double from z0 = 150 m (ASTRA) / 250 m (TOAD). | safety | **Fixed by FIX-F07** (shift-invariant reparametrisation; max cond 2.8e3). |
| F08 | `Guidance/TrajectoryManager.m:43`–`:67` | `AbortGuidance.cpp:67` | `Report.Feasible` and PolyQP errors were ignored; there was no fallback. | safety | **Fixed by FIX-F08** (fallback descent, `AbortStatus`; D1, D2). The C++ "non-finite estimate → fallback" branch was removed (NO-NAN, Phase 3); MATLAB keeps it. |
| F09 | `Guidance/TrajectoryManager.m:14` | `AbortGuidance.cpp:10` | The trigger is a hard-coded constant (`ABORT = 0`, never fires). See F31. | safety | **Moot (ASTRA integration Phase 6)**: the trigger is removed; the caller decides when to abort and calls `Abort::plan` once. |
| F10 | `Guidance/PolynomialQP/PolyQP.m:3`, `:25` | `PolyQP.cpp:13`, `:38` | The `nargin < 2` / `~Limits.HasVehicle` branches are dead. Transpiled anyway. | info | Open (the C++ has no `Constants` argument: the ASTRA constants are compiled in). |
| F11 | `Guidance/PolynomialQP/evaluateConstraints.m:6`, `sampleOutputs.m:5` | `evaluateConstraints.cpp:16`, `sampleOutputs.cpp:11` | The 3999 caps silently truncate the grids for long trajectories. | behaviour | Open. For ASTRA the output grid is capped lower, at `POLYQP_MAX_NODES` = 601 (a longer plan returns `CapacityExceeded` and flies the fallback); the evaluation grid is streamed and keeps the 3999 cap. |
| F12 | `Guidance/PolynomialQP/evaluateConstraints.m:128` | `evaluateConstraints.cpp:150` | Dilation exponents differ from paper §3.5.2 (T⁻³). The code is followed. | info | Open. |
| F13 | `Guidance/PolynomialQP/PolyEval.m:21` | `PolyEval.cpp:23` | `hasCruise` is always true for a `QPTraj` from `solveAllLegs`; transpiled as constant `true`. | info | Open. |
| F14 | `Guidance/PolynomialQP/assembleQP.m:2`, `solveLeg.m:2` | `assembleQP.cpp:9`, `solveLeg.cpp:9` | "Paper Eq. ??": equation numbers are missing in the source comments. Preserved verbatim. | info | Open. |
| F15 | `Helper/DCM_Quat_Conversion.m` | `cpp/matlab_funcs.h` (double overload) | Duplicates the `float` version in `cpp/matlab_helpers.cpp`; the two must stay in sync. | info | Open. Now a double overload beside the float one (Phase 1); `DCM_Quat_Conversion.cpp` is deleted. |
| F16 | `Guidance/TrajectoryManager.m:48` | `AbortGuidance.cpp:68` | `planStatus` and `Report` are returned to the C++ caller. | info | Open (by design); `Abort::Result::solver_status` carries the status to the flight log. |

## Flags found during the transpile (F24–F33)

| ID | File:line (.m) | File:line (.cpp) | Description | Severity | Status |
| --- | --- | --- | --- | --- | --- |
| F24 | `resolveWaypoints.m:3-6`, `solveLeg.m:8`, `sampleOutputs.m:6` | `PolyQP_config.h:28` | Derived capacities: `POLYQP_MAX_NODES = POLYQP_MAX_GRID + 1`, `POLYQP_MAX_NV = POLYQP_MAX_R * POLYQP_MAX_KNOTS`. | info | Open. |
| F25 | n/a | `PolyQP_config.h:37` | Fixed-max Eigen buffers exceed Eigen's 128 KB static-size assert, so `EIGEN_STACK_ALLOCATION_LIMIT` is 0. Every such buffer is `static`. | info | **Removed (Phase 4)**: at the abort capacities every buffer is small and fully fixed-size, so the flight build leaves Eigen's limit alone (only the verification MEX `polyqp_mex_wide` sets it). |
| F26 | n/a | `matlab_compat.h:4` | The fast-math guard tests the value of `__FINITE_MATH_ONLY__`, not only `defined()`. | info | **Removed (Phase 3)** with the guard and `matlab_compat.h`. The solver no longer uses NaN, but it is verified only without fast math: build with `-ffp-contract=off` and no `-ffast-math` (`cpp/INTEGRATION.md`). |
| F27 | `Guidance/PolynomialQP/flatnessInvert.m:54` | `flatnessInvert.cpp:75` | The `x_C_alt` fallback changed b2, but `cd`/`cdd` still differentiated `cross(u, x_C)`. | behaviour | **Fixed by FIX-F27**. |
| F28 | `Guidance/TrajectoryManager.m:27` | `AbortGuidance.cpp:21` | The comment says "Minimum-jounce" (r = 4) but `R = 5` is minimum-crackle. | info | **Removed from the C++ (Phase 2)** with `Opts`: r = 5 is the constexpr `POLYQP_R`. The MATLAB comment is unchanged. |
| F29 | `Guidance/PolynomialQP/flatnessInvert.m:109` | `flatnessInvert.cpp:136` | Propellant integration started from the full-tank masses, not `X_est(14:15)`. | behaviour | **Removed (vehicle 1, Phase 2)**: FIX-F29 is superseded, ASTRA carries no propellant. The C++ has no propellant integration; `X(14:15)` stays 0. |
| F30 | `Guidance/PolynomialQP/allocateControls.m:24`, `:36` | `allocateControls.cpp:54` | U(4) is a roll torque in N·m. That is consistent with the plant and the LQR design model, but the controller labels it rad/s² and clamps it at ±4. `evaluateConstraints` compares it with `RollLimit = 0.9 * MaxRollRate` (4.5 N·m). | safety | **BLOCKER B1** (controller side; reported, not edited). |
| F31 | `Guidance/TrajectoryManager.m:14-21` | `AbortGuidance.cpp:11`, `:12` | `ABORT` is a parameter defaulting to 0, so flight code cannot abort unless a caller deliberately sets it; tests set it to 1. The persistents `tAbort`/`abortCalculated`/`abortStatus` live in `AbortGuidanceState`. | info | **Moot (Phase 6)**: `ABORT`, `tAbort`, `abortCalculated` and `AbortGuidanceState` are removed in MATLAB and C++; the caller decides when. |
| F32 | `Guidance/TrajectoryManager.m:74` | `AbortGuidance.cpp:69` | On a PolyQP error the C++ used to return before `abortCalculated` was set, so the abort was recomputed every tick. | behaviour | **Moot (Phase 6)**: the one-attempt latch went with the trigger; each `Abort::plan` call plans exactly once. |
| F33 | `Guidance/PolynomialQP/PolyQP.m:6-13` | `PolyQP.cpp:14` | `nargin < 3 / isempty(Opts)` is the default argument `Opts = PolyQPOpts()`. | info | Open; `Opts`/`PolyQPOpts` are removed from the C++ (r = 5 hardcoded, Phase 2). |

## Flags added on `abort-fixes` (F34–F38)

| ID | File:line (.m) | File:line (.cpp) | Description | Severity |
| --- | --- | --- | --- | --- |
| F34 | `Guidance/PolynomialQP/abortFallback.m:35` | `abortFallback.cpp:41` | The fallback's 3999-sample cap is also clamped to the sink capacity (at most `POLYQP_MAX_NODES`), so the fallback always fits the caller's arrays; a long descent is sampled more coarsely. Identical to MATLAB when it fits. | info |
| F35 | `Guidance/TrajectoryManager.m:46` | `AbortGuidance.cpp:75` | The C++ catches every PolyQP error, including `CapacityExceeded`, through the status return. MATLAB catches PolyQP errors only in interpreted MATLAB (`coder.target('MATLAB')`): code generation has no `try/catch`, so in Simulink a PolyQP error still stops the simulation. | safety |
| F36 | `Helper/quatRot.m` | `cpp/matlab_funcs.h` (double overload) | Duplicates the `float` `quatRot` in `cpp/matlab_helpers.cpp` in double; the two must stay in sync. Added for FIX-F03; `quatRot.cpp` is deleted (Phase 1). | info |
| F37 | `Guidance/TrajectoryManager.m:16` | `AbortGuidance.cpp:13` | `ABORT_PIN_THRUST` is a parameter defaulting to the MATLAB constant (false), like `ABORT` (F31), so tests can switch it on. | info — **removed (Phase 5)**: the start BC is always on. |
| F38 | `Guidance/PolynomialQP/buildVehicleLimits.m:1`, `:21` | `buildVehicleLimits.cpp:9` | The optional MATLAB `Opts` (`nargin`, `isfield`, `isempty`) becomes the default argument `Opts = PolyQPOpts()`, placed last. A missing, empty or non-finite `InitOxMass`/`InitFuMass` is the NaN default (full tank). | info — **superseded (Phase 2)**: ASTRA carries no propellant; `Opts` and the `Constants` argument are not carried. |

## Deviations from the MATLAB

All are permitted by spec Section 5.8 (index conversion, `error` to status return, out-parameters, fixed-max storage with capacity checks, `static` large buffers, `arguments` to defaults, `char[10]` labels, field-by-field `struct` literals, float↔double at the `Abort::plan` boundary), except those flagged above: F24, F32 (fixed, now moot), F33, F34, F35, F38. Numerically there are no deviations: the linear algebra follows spec Section 5.2. `findFirst` (`cpp/matlab_funcs.h`, FIX-F07) is MATLAB's `find(idx == v)`.

The ASTRA integration added these, each tagged in the C++ (the self-review in `verify_polyqp` checks that every MATLAB comment still sits on its code):

- `ASTRA-ONLY`: vehicle 1 and r ≠ 5 removed; the constants are `constexpr` in `ASTRA_PolyQP_Constants.h`; `ComputeJtot` reduces to `J_tot = J`, `CGz = rTB`.
- `NO-NAN`: a MATLAB NaN "free" boundary condition is a parallel `bool` "fixed" mask (free entries hold 0); `isnan` tests read the mask. No NaN is created anywhere in `PolyQP/`; the MEX converts NaN ↔ masks at the verification boundary.
- `STATIC`: every Eigen type is fully fixed-size (`cpp/matlab_funcs.h`); `R_PP` is solved as a `POLYQP_MAX_FREE` × `POLYQP_MAX_FREE` system padded with an identity block, with the rcond guard estimated on the unpadded block.
- `STREAM`: `evaluateConstraints` streams over its grid (per-point `flatnessInvertPoint`, previous `U` and quaternion carried); `sampleOutputs` and `abortFallback` write each output node straight into the caller's arrays through a `NodeSink`.

## Non-reentrancy (static buffers)

These buffers are `static`, so **no PolyQP function, and not `Abort::plan`, is reentrant or thread-safe**. Calls must be serialised; do not call them from an ISR or from two tasks.

| File | Static buffers |
| --- | --- |
| `PolyQP.cpp` | `g_workspace` (`PolyQPWorkspace`, 24,880 B), shared by every stage, which run sequentially: `solveLeg` (`Rmat`, `Hseg`, `Rns`, `fx`/`fxFixed`, `dv`, `dF`, `R_PP`, `R_PF`, `Tz`, `Rns_PP`, `TzT_Rns_PP`, `R_PP_T`, `Rns_PF`, `d_diag`, `s_jac`, `s_jac_s_jacT`, `R_PP_eq`, `R_PP_eq_T`, `R_PF_dF`, `rhs_eq`, the `LLT`/`PartialPivLU` workspaces, `R_PP_eq_rhs_eq`, `z`, `Tz_z`); `solveAllLegs` (`FixP`, `CoefAll`, `DerivAll`, `FixLegAT`, `FixLeg`, `CoefLeg`, `DerivLeg`, `FixOut`, with masks); `allocateTimes` (`CoefIgnored`, `DerivIgnored`, `Deriv`); `resolveWaypoints` (`FixCg`); `PolyQP` (`Plan`, `EvalOut`, `dLeg`); `AbortGuidance` (`Waypoints`, `QPTraj`) |
| `assembleQP.cpp` | 32 B of function-local constants |
| `../Abort.cpp` | `report` (`ReportT`, 144 B) |
| Verification only (`POLYQP_ERROR_MESSAGES`) | `lastError()`, `lastRcond()` in `PolyQP_types.h` |

The grid functions `flatnessInvert`, `allocateControls` and `PolyEval` keep their MATLAB outputs for verification; the flight path does not call them on full grids.

## Memory footprint (measured, host build, `size -m`, `.bss` + `.data`)

| Profile | `abort-cpp` | `abort-fixes` | `abort-astra-integration` |
| --- | --- | --- | --- |
| MATLAB capacities (`POLYQP_MAX_SEG=40`, `POLYQP_MAX_GRID=3999`) | 9.58 MB | 13.04 MB | n/a (ASTRA only; the verification MEX `polyqp_mex_wide` builds the same sources with 16 segments and a 3999 grid) |
| Abort-sized | 1.67 MB (`-DPOLYQP_MAX_SEG=4 -DPOLYQP_MAX_GRID=1200`) | 1.81 MB (same) | **25,056 B** (defaults: 3 waypoints, 3 segments, 4 knots, 5 free, 601 output nodes) |

`abort-astra-integration`, per object (flight flags, `build_host.sh` step 2b): `PolyQP.o` 24,880 B, `Abort.o` 144 B, `assembleQP.o` 32 B, every other object 0 B. The target is ≤ 256 KB. The caller's abort trajectory buffer (`Abort::MAX_POINTS` × 80 B = 48,080 B) is outside PolyQP; see `cpp/INTEGRATION.md`. Host stack at the deepest call chain is about 13 KB (`-fstack-usage`).
