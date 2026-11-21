# PolyQP transpile: flag log

Every `TRANSPILE-FLAG` in `cpp/PolyQP/` is listed here. A flag marks something in the MATLAB that looks wrong, or a place where the C++ had to differ from the MATLAB beyond the permitted deviations (spec Section 5.8).

On branch `abort-cpp` nothing flagged was fixed: the C++ did what the MATLAB at `5ac6bd23` does. On branch `abort-fixes` the safety-relevant flags were fixed **in the MATLAB first and mirrored in the C++** (see `ABORT_FIXES.md`). A fixed flag's comment reads `TRANSPILE-FLAG Fxx (fixed by FIX-Fxx)`. Every changed line carries a `FIX-Fxx` comment in both languages. The **Status** column below gives the state on `abort-fixes`.

MATLAB paths are relative to `TOAD/`; C++ paths to `TOAD/cpp/PolyQP/`. Line numbers are at `abort-fixes`. IDs F17–F23 are unused: new flags were numbered from F24.

Severity: **info** (no effect on outputs), **behaviour** (changes what the planner produces compared with what was presumably intended), **safety** (can affect the vehicle in an abort).

## Pre-identified flags (F01–F16)

| ID | File:line (.m) | File:line (.cpp) | Description | Severity | Status |
| --- | --- | --- | --- | --- | --- |
| F01 | `Guidance/TrajectoryManager.m:15`, `:41` | `AbortGuidance.cpp:14`, `:48` | Abort waypoint 2 altitude is 2 m; the lead's spec says 1 m above the pad. | behaviour | **Open, DECISION-D5**: now the named constant `ABORT_HOVER_ALT = 2` (no value change). |
| F02 | `Guidance/TrajectoryManager.m:16`, `:37` | `AbortGuidance.cpp:42` | By default only position and velocity are pinned at knot 1. Accel/jerk/snap are free, so thrust and attitude feedforward can step at the abort instant (measured: up to 25° attitude, 41 % thrust). | safety | **Option added (FIX-F02), off (DECISION-D4)**: `ABORT_PIN_THRUST` pins accel/jerk to hover along the current b3. It removes the step, but makes tilted aborts infeasible (lateral drift). |
| F03 | `Guidance/PolynomialQP/PolyQP.m:140`, `:223` | `PolyQP.cpp:169`, `:270` | Knot-1 P-offset assumed b3 = e3 and was never refined. The velocity BC was applied to dP/dt, not v_cg + rho·db3/dt. | behaviour | **Fixed by FIX-F03** (`B3BC`/`B3DotBC`, abort passes the estimated attitude). |
| F04 | `Guidance/PolynomialQP/PolyQP.m:186` | `PolyQP.cpp:228` | `nPasses = 3`: for the abort, passes 2–3 recompute pass 1 exactly. Wasted work in flight. | info | Open (out of scope). |
| F05 | `Guidance/PolynomialQP/PolyQP.m:37` | `PolyQP.cpp:62` | `legCanContract` is false for every abort, so the abort never contracts below its initial time guess. | behaviour | Open (out of scope). |
| F06 | `Guidance/PolynomialQP/resolveWaypoints.m:90`, `:96` | `resolveWaypoints.cpp:108`, `:116` | `all(~isnan(dp))` is false with NaN lateral positions, so `StraightLine` on abort waypoint 3 is inert and the abort never takes the cruise split. | behaviour | **Consequence addressed by FIX-F07**: the free lateral landing position no longer degrades conditioning. `StraightLine` itself is still inert. |
| F07 | `Guidance/PolynomialQP/solveLeg.m:94` | `solveLeg.cpp:138` | Lateral R_PP was near-singular for unbalanced segment times (cond ≈ (T1/T2)^9). `SingularRPP` tripped in double from z0 = 150 m (ASTRA) / 250 m (TOAD). | safety | **Fixed by FIX-F07** (shift-invariant reparametrisation; max cond 2.8e3). |
| F08 | `Guidance/TrajectoryManager.m:43`–`:67` | `AbortGuidance.cpp:67` | `Report.Feasible` and PolyQP errors were ignored; there was no fallback. | safety | **Fixed by FIX-F08** (fallback descent, `AbortStatus`; D1, D2). |
| F09 | `Guidance/TrajectoryManager.m:14` | `AbortGuidance.cpp:10` | The trigger is a hard-coded constant (`ABORT = 0`, never fires). See F31. | safety | Open (deliberate; flight trigger is an integration decision). |
| F10 | `Guidance/PolynomialQP/PolyQP.m:3`, `:25` | `PolyQP.cpp:13`, `:38` | The `nargin < 2` / `~Limits.HasVehicle` branches are dead. Transpiled anyway. | info | Open. |
| F11 | `Guidance/PolynomialQP/evaluateConstraints.m:6`, `sampleOutputs.m:5` | `evaluateConstraints.cpp:16`, `sampleOutputs.cpp:11` | The 3999 caps silently truncate the grids for long trajectories. Depleted-tank TOAD aborts now reach 400+ s (FIX-F29), i.e. a coarser evaluation grid. | behaviour | Open. |
| F12 | `Guidance/PolynomialQP/evaluateConstraints.m:128` | `evaluateConstraints.cpp:150` | Dilation exponents differ from paper §3.5.2 (T⁻³). The code is followed. | info | Open. |
| F13 | `Guidance/PolynomialQP/PolyEval.m:21` | `PolyEval.cpp:23` | `hasCruise` is always true for a `QPTraj` from `solveAllLegs`; transpiled as constant `true`. | info | Open. |
| F14 | `Guidance/PolynomialQP/assembleQP.m:2`, `solveLeg.m:2` | `assembleQP.cpp:9`, `solveLeg.cpp:9` | "Paper Eq. ??": equation numbers are missing in the source comments. Preserved verbatim. | info | Open. |
| F15 | `Helper/DCM_Quat_Conversion.m` | `DCM_Quat_Conversion.cpp:2` | Duplicates the `float` version in `cpp/matlab_helpers.cpp`; the two must stay in sync. | info | Open. |
| F16 | `Guidance/TrajectoryManager.m:48` | `AbortGuidance.cpp:68` | `state.planStatus` and `Report` are returned to the C++ caller. MATLAB keeps only the persistent `abortStatus` (the Simulink block's outputs are unchanged). | info | Open (by design). |

## Flags found during the transpile (F24–F33)

| ID | File:line (.m) | File:line (.cpp) | Description | Severity | Status |
| --- | --- | --- | --- | --- | --- |
| F24 | `resolveWaypoints.m:3-6`, `solveLeg.m:8`, `sampleOutputs.m:6` | `PolyQP_config.h:28` | Derived capacities: `POLYQP_MAX_NODES = POLYQP_MAX_GRID + 1`, `POLYQP_MAX_NV = POLYQP_MAX_R * POLYQP_MAX_KNOTS`. | info | Open. |
| F25 | n/a | `PolyQP_config.h:37` | Fixed-max Eigen buffers exceed Eigen's 128 KB static-size assert, so `EIGEN_STACK_ALLOCATION_LIMIT` is 0. Every such buffer is `static`. | info | Open. |
| F26 | n/a | `matlab_compat.h:4` | The fast-math guard tests the value of `__FINITE_MATH_ONLY__`, not only `defined()`. | info | Open. |
| F27 | `Guidance/PolynomialQP/flatnessInvert.m:54` | `flatnessInvert.cpp:75` | The `x_C_alt` fallback changed b2, but `cd`/`cdd` still differentiated `cross(u, x_C)`. | behaviour | **Fixed by FIX-F27**. |
| F28 | `Guidance/TrajectoryManager.m:27` | `AbortGuidance.cpp:21` | The comment says "Minimum-jounce" (r = 4) but `R = 5` is minimum-crackle. | info | Open. |
| F29 | `Guidance/PolynomialQP/flatnessInvert.m:109` | `flatnessInvert.cpp:136` | Propellant integration started from the full-tank masses, not `X_est(14:15)`. | behaviour | **Fixed by FIX-F29** (`Opts.InitOxMass`/`InitFuMass` through `buildVehicleLimits`). |
| F30 | `Guidance/PolynomialQP/allocateControls.m:24`, `:36` | `allocateControls.cpp:54` | U(4) is a roll torque in N·m. That is consistent with the plant and the LQR design model, but the controller labels it rad/s² and clamps it at ±4. `evaluateConstraints` compares it with `RollLimit = 0.9 * MaxRollRate` (4.5 N·m). | safety | **BLOCKER B1** (controller side; reported, not edited). |
| F31 | `Guidance/TrajectoryManager.m:14-21` | `AbortGuidance.cpp:11`, `:12` | `ABORT` is a parameter defaulting to 0, so flight code cannot abort unless a caller deliberately sets it; tests set it to 1. The persistents `tAbort`/`abortCalculated`/`abortStatus` live in `AbortGuidanceState`. | info | Open (by design). |
| F32 | `Guidance/TrajectoryManager.m:74` | `AbortGuidance.cpp:69` | On a PolyQP error the C++ used to return before `abortCalculated` was set, so the abort was recomputed every tick. | behaviour | **Fixed by FIX-F32**. |
| F33 | `Guidance/PolynomialQP/PolyQP.m:6-13` | `PolyQP.cpp:14` | `nargin < 3 / isempty(Opts)` is the default argument `Opts = PolyQPOpts()`. | info | Open. |

## Flags added on `abort-fixes` (F34–F38)

| ID | File:line (.m) | File:line (.cpp) | Description | Severity |
| --- | --- | --- | --- | --- |
| F34 | `Guidance/PolynomialQP/abortFallback.m:35` | `abortFallback.cpp:41` | The fallback's 3999-sample cap is also clamped to `POLYQP_MAX_GRID`, so the fallback cannot exceed a reduced flight capacity profile. Identical to MATLAB at the default capacities. | info |
| F35 | `Guidance/TrajectoryManager.m:46` | `AbortGuidance.cpp:75` | The C++ catches every PolyQP error, including `CapacityExceeded`, through the status return. MATLAB catches PolyQP errors only in interpreted MATLAB (`coder.target('MATLAB')`): code generation has no `try/catch`, so in Simulink a PolyQP error still stops the simulation. | safety |
| F36 | `Helper/quatRot.m` | `quatRot.cpp:8` | Duplicates the `float` `quatRot` in `cpp/matlab_helpers.cpp` in `Real`; the two must stay in sync. Added for FIX-F03. | info |
| F37 | `Guidance/TrajectoryManager.m:16` | `AbortGuidance.cpp:13` | `ABORT_PIN_THRUST` is a parameter defaulting to the MATLAB constant (false), like `ABORT` (F31), so tests can switch it on. | info |
| F38 | `Guidance/PolynomialQP/buildVehicleLimits.m:1`, `:21` | `buildVehicleLimits.cpp:9` | The optional MATLAB `Opts` (`nargin`, `isfield`, `isempty`) becomes the default argument `Opts = PolyQPOpts()`, placed last. A missing, empty or non-finite `InitOxMass`/`InitFuMass` is the NaN default (full tank). | info |

## Deviations from the MATLAB

All are permitted by spec Section 5.8 (index conversion, `error` to status return, out-parameters, `persistent` to a state struct, fixed-max storage with capacity checks, `static` large buffers, `arguments` to defaults, `char[10]` labels, field-by-field `struct` literals, float↔double at the `AbortGuidance` boundary), except those flagged above: F24, F25, F26, F31, F32 (fixed), F33, F34, F35, F37, F38. Numerically there are no deviations: the linear algebra follows spec Section 5.2. `findFirst` (`matlab_compat.h`, FIX-F07) is MATLAB's `find(idx == v)`.

## Non-reentrancy (static buffers)

These buffers are `static` (deviation 6), so **no PolyQP function is reentrant or thread-safe**. Calls must be serialised; do not call PolyQP from an ISR or from two tasks.

| File | Static buffers |
| --- | --- |
| `solveLeg.cpp` | `Rmat`, `R_PP`, `R_PF`, `s_jac_s_jacT`, `R_PP_eq`, `R_PP_eq_T`, `Hseg`; FIX-F07: `Tz`, `Rns`, `Rns_PP`, `TzT_Rns_PP`, `R_PP_T`, `Rns_PF` (nv × nv each), the `PartialPivLU`/`LLT` workspaces `lu_rc`, `llt`, `lu`, and the nv-vectors `fx`, `dv`, `dF`, `d_diag`, `s_jac`, `R_PF_dF`, `rhs_eq`, `R_PP_eq_rhs_eq`, `z`, `Tz_z` |
| `flatnessInvert.cpp` | `P_pos`, `P_vel`, `P_acc`, `P_jerk`, `P_snap`, `q_s`, `b1_s`, `b2_s`, `b3_s`, `b3d_s`, `b3dd_s`, `w_s`, `wd_s`, `n_s`, `r_s`, `v_s`, `ml_s`, `mi_s`, `F_body`, `Tau_body`, `CGz_vec` |
| `evaluateConstraints.cpp` | `t_grid`, `X`, `Wrench`, `AccelCg`, `U`, `U_rate`, `diffU`, `seg_idx` |
| `sampleOutputs.cpp` | `t_vec`, `X_nodes`, `Wrench`, `AccelCgIgnored`, `U_nodes`, `U_ctrl` |
| `PolyEval.cpp` | `idx`, `t_loc`, `tp`, `c_row` |
| `PolyQP.cpp` | `Plan`, `EvalOut`, `dLeg`; in `solveAllLegs`: `FixP`, `CoefAll`, `DerivAll`, `midVelBins`, `FixLeg` (×2), `CoefLeg`, `DerivLeg`, `FixOut` |
| `allocateTimes.cpp` | `CoefIgnored`, `DerivIgnored`, `Deriv` |
| `resolveWaypoints.cpp` | `FixCg` |
| `AbortGuidance.cpp` | `Waypoints`, `InitialGuess`, `QPTraj`; float overload: `Time_real`, `States_real`, `Inputs_real` |
| `abortFallback.cpp` | `t_vec` |

## Memory footprint (measured, host build, `size -m`, bss)

| Profile | `abort-cpp` | `abort-fixes` |
| --- | --- | --- |
| MATLAB capacities (`POLYQP_MAX_SEG=40`, `POLYQP_MAX_GRID=3999`) | 9.58 MB | **13.04 MB** (`solveLeg` 7.33 MB, `flatnessInvert` 1.79 MB, `AbortGuidance` 1.33 MB, `evaluateConstraints` 1.23 MB, `sampleOutputs` 1.09 MB) |
| Abort-sized (`-DPOLYQP_MAX_SEG=4 -DPOLYQP_MAX_GRID=1200`) | 1.67 MB | **1.81 MB** |

The growth is FIX-F07's six extra nv × nv buffers in `solveLeg` (nv = 246 at MATLAB capacities, 25 at the abort profile). Several could share storage with existing buffers; the memory profile was out of scope for the fixes, so this is reported, not optimised. Neither profile fits the ~1 MB of on-chip RAM on a Teensy 4.x or STM32H7, so a flight build still needs a smaller capacity profile and/or external RAM. An abort that exceeds a reduced profile now returns `CapacityExceeded` inside `AbortGuidance`, which flies the fallback (`AbortStatus::FALLBACK`, F35).
