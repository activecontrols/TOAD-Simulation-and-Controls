# Abort fixes (branch `abort-fixes`)

These fixes address the problems on branch `abort-cpp` that would break an abort in flight. The MATLAB is still the specification: each fix was made in the `.m` first and then mirrored in the C++ with the transpile's 1:1 conventions. Every changed line carries `% FIX-Fxx: <reason>` in MATLAB and the identical `// FIX-Fxx: <reason>` in C++; the automated comment check in `verify_polyqp.m` enforces the correspondence. Fixed flags read `TRANSPILE-FLAG Fxx (fixed by FIX-Fxx)`; `TRANSPILE_FLAGS.md` has the status of every flag.

Paths: MATLAB relative to `TOAD/`, C++ relative to `TOAD/cpp/PolyQP/`. Line numbers are at `abort-fixes`. All numbers below were measured by the scripts listed under "Reproducing". Nothing here is extrapolated.

**Summary.** Six fixes were made in MATLAB and C++: F07, F08/F32, F01 (D3), F03, F29 and F27. One option was added and left off: F02 (D4). There is one blocker for the controller owners, B1 (roll units). Every MATLAB/C++ parity check passes at every tier, all seven regression sections pass, the abort block still generates code, and the host smoke test and flight-flags `-Werror` build pass.

## Before / after

Case sets:
- **Frozen:** the 200 + 200 random abort states of `verify_polyqp.m` (seed 42), upright, zero rates, full tanks.
- **Sweep:** start at (5, −3, z0) with z0 from 0.5 to 300 m and 4 velocities within the caps (`diagnose_aborts.m` §1.2), plus 300 further runs in `fix_regressions` F07.
- **Tilted:** 100 per vehicle, tilt ≤ 30°, |ω| ≤ 1 rad/s.

| Metric | `abort-cpp` (baseline) | `abort-fixes` |
| --- | --- | --- |
| `SingularRPP` in the conditioning sweep | ASTRA 32 / 116 (every run from z0 = 150 m), TOAD 12 / 116 (from 250 m) | **0 / 116** each; 0 / 300 in the F07 regression sweep |
| Max cond(R_PP_eq) over successful solves | 7.6e12 (ASTRA), 8.0e12 (TOAD); 6.7e15 at the initial times at 300 m | **2.8e3**, independent of altitude (fit exponent on T1/T2: 9.0 before, 0.0 after) |
| Infeasible, frozen ASTRA | 82 / 200 (41.0 %): all 82 over-cap start velocity | **0 / 200** |
| Infeasible, frozen TOAD | 93 / 200 (46.5 %): all 93 over-cap start velocity | **0 / 200** |
| Infeasible, within-cap / over-cap / tilted sets (100 each per vehicle) | – (new sets) | **0 / 0 / 0** for both vehicles |
| Infeasible, within-cap starts (frozen) | 0 | 0 |
| P-space start BC error (tilted starts) | ρ\|b3 − e3\|: up to 0.11 m (ASTRA), 0.31 m (TOAD); velocity ρ\|ω × b3\| up to 0.22 / 0.60 m/s | **0** (exact) |
| Plan CG at t = 0 vs `X_est`, upright frozen, max | 0.097 m / 0.135 m/s (ASTRA), 0.181 m / 0.068 m/s (TOAD) | 0.090 m / 0.135 m/s, 0.175 m / 0.068 m/s (pin off); ≤ 4e-15 (pin on) |
| Plan CG at t = 0 vs `X_est`, tilted, max (pin off) | 0.092 m / 0.064 m/s (ASTRA), 0.145 m / 0.132 m/s (TOAD) | 0.164 m / 0.204 m/s, 0.336 m / 0.561 m/s (see D6); ≤ 2e-15 with the pin |
| Attitude step at the abort instant, upright frozen, max | 25.5° (ASTRA), 17.3° (TOAD) | 23.6°, 16.7° (pin off, default); 0° (pin on) |
| PolyQP error in flight | C++ returned with no trajectory and retried every tick | One attempt; closed-form fallback (`AbortStatus::FALLBACK`) |
| Float build, 200 ASTRA aborts | 10 status failures; position error max 26.6 m; T_total error max 130 s | **0** failures; position max 0.021 m; T_total max 0.018 s |
| Float `solveLeg`, abort leg at (10, 0.6) s and (8, 0.8) s | `SingularRPP` | solves; Deriv rel err 8.1e-6 and 2.0e-6 |
| C++ vs MATLAB, `AbortGuidance` end-to-end (abs) | 2.2e-5 max (21 strict misses, round-off-explained) | 2.9e-9 max (0 strict misses) |
| Static RAM (bss), MATLAB capacities / abort profile | 9.58 MB / 1.67 MB | 13.04 MB / 1.81 MB (FIX-F07 buffers; memory was out of scope) |

## Fixes

### FIX-F07 (covers F06's consequence): shift-invariant segments in `solveLeg`

- **Cause.** The abort's landing x, y appear twice as free unknowns: p2 at ABT GND and p3 at ABT END. Segment 2 has every non-position derivative fixed at both ends, so its cost depends only on p3 − p2. The direction p2 + p3 is fixed only by the long segment 1, and cond(R_PP_eq) grows as (T1/T2)^9. Measured: an exponent of 8.99, with the guard predicted to trip at T1/T2 = 14.2. It actually tripped at z0 = 150 m on ASTRA (T1/T2 = 14.8; 140 m still passed) and 250 m on TOAD, whose T2 is 7.9 s against ASTRA's 4.9 s. The earlier "≈ 140 m" estimate holds for ASTRA only.
- **Change.** For each axis and segment, if derivatives 1..r−1 are fixed at both knots and both positions are free, the free unknown p_{m+1} is replaced by Δ_m = p_{m+1} − p_m (d_P = Tz·z, which handles chains). The existing Jacobi equilibration and rcond guard are applied to the reduced matrix, and every other case is untouched.
  - `Guidance/PolynomialQP/solveLeg.m:13`, `:20`, `:38-80`, `:100-102`, `:106-117`
  - `solveLeg.cpp:23`, `:36`, `:55-119`, `:153-159`, `:163-181`
  - `matlab_compat.h:148` (`findFirst`, MATLAB `find(idx == v)`)
- **Deviation from the prompt's formula (DECISION D8).** Computing Tzᵀ R_PP Tz from the already summed R_PP is exact in exact arithmetic but not in floating point. The eps·T2⁻⁹ rounding in the large segment-2 entries survives the cancellation. Measured: a 1.4e-1 relative error against `vpa` at T = (9.96, 0.27) s, a case the old code at least rejected with `SingularRPP`, and float still failed at (10, 0.6) s. The reduced system is therefore assembled segment by segment: a shift-invariant segment contributes only its Δ_m row and column (local index r + 1), and its p_m part is dropped. That part is identically zero in exact arithmetic. The cost is evaluated the same way. The minimiser is the same.
- **Acceptance** (`fix_regressions` F07):
  - **(a)** 1,250 baseline solves that succeeded: non-shift cases are bit-identical in MATLAB. Of the shift cases, 252 differ from the baseline by more than 1e-9. In every one, the fixed result is at least as close to the 60-digit `vpa` solve as the baseline was, so the difference is the baseline's own rounding. Error against `vpa` on those cases, median / 95th / max:
    - baseline: 6.0e-9 / 9.3e-6 / 1.6e-4;
    - fixed MATLAB: 6.1e-14 / 1.8e-7 / 1.6e-4;
    - fixed C++: 6.7e-14 / 1.3e-7 / 8.5e-5.

    The maxima are non-shift mission legs at cond ≈ 1e9, which are unchanged.
  - **(b)** 0 `SingularRPP` to z0 = 300 m. Max cond(R_PP_eq) is 2.77e3 (expected ≈ 3e3), independent of altitude.
  - **(c)** 204 of 221 baseline error cases keep the same ID. The other 17 now solve: each had a shift-invariant segment and was `SingularRPP` only from the degenerate direction, and each now matches `vpa` to ≤ 1.4e-11. No ID changed otherwise.
  - **(d)** Float: 0 status failures (was 10); the float numbers are in the table above.
  - Code generation found `Rns` undefined on some paths; fixed in `e9c234f`.

### FIX-F08 / FIX-F32: defined behaviour when the abort plan fails

- **Cause.** MATLAB ignored `Report.Feasible` and PolyQP errors. On an error the C++ returned before `abortCalculated` was set, so it recomputed every tick and produced no trajectory.
- **Change.**
  - `abortCalculated` is set on every path (`Guidance/TrajectoryManager.m:74`, `AbortGuidance.cpp:98`).
  - A non-finite `X_est(5:10)`, or a PolyQP error, gives `abortFallback` (D1):
    - `Guidance/PolynomialQP/abortFallback.m` (new), `abortFallback.cpp` (new, flag F34);
    - `Guidance/TrajectoryManager.m:43-70`, `AbortGuidance.cpp:65-95`.
  - An infeasible plan is flown and flagged (D2).
  - MATLAB keeps a persistent `abortStatus`: 0 OK, 1 INFEASIBLE_FLOWN, 2 FALLBACK (`TrajectoryManager.m:19`, `:22`). The C++ returns `AbortStatus` and keeps `state.abortStatus`, plus `state.planStatus` (PolyQP's status) for logging (`AbortGuidance.h:18-38`). The `TrajectoryManager` outputs are unchanged, so the Simulink block is unaffected.
- **Code-generation constraint (DECISION D7).** `TrajectoryManager.m` runs in a Simulink MATLAB Function block (`Guidance/TOAD_TrajectoryGen_3.slx`), and MATLAB code generation has no `try/catch` (measured on R2026b). The `try/catch` around `PolyQP` is therefore inside `if coder.target('MATLAB')`. Interpreted MATLAB (the verification harness) falls back on any PolyQP error. In Simulink, a PolyQP error still stops the simulation (F35); only the non-finite-input pre-check routes to the fallback there. The flight C++ catches every error, including `CapacityExceeded`, through the status return. `verify/check_abort_codegen.m` confirms the block generates code.
- **Acceptance** (`fix_regressions` F08):
  - A forced PolyQP error (z0 = 1e300 → `SingularRPP`) and a NaN velocity both produce the fallback in MATLAB and C++. A second call recomputes nothing.
  - On all 400 frozen states plus 4 edge cases (on the ground, below ground, 300 m, non-finite position and propellant) the fallback starts exactly at the estimated position and ends at rest at z = 0 exactly. |vz| stays ≤ 1 m/s, samples are within cap × dt, and the output is upright, hover thrust and finite.
  - MATLAB/C++ parity holds on all three paths for both vehicles (planned paths ≤ 1.5e-9 abs; fallback within 1.4 × MATLAB's own `vpa`-measured rounding).

### FIX-F01 (DECISION D3): over-cap start velocity

- **Cause (Phase 1.1, hypothesis confirmed).** Segment 1's velocity cap is waypoint 2's default `MaxVel = [4; 4; 5]` m/s. A start velocity pinned above it is infeasible from t = 0, and time dilation cannot reduce a pinned velocity, so the outer loop stalls.
  - Predicted over-cap rate 42.4 %; measured 42.5 % (ASTRA) and 47.5 % (TOAD).
  - Every infeasible case was over-cap: 82 / 82 and 93 / 93.
  - In 173 of 175 velocity dominated, in segment 1, from t = 0. In the other 2 (ASTRA) acceleration dominated, also from t = 0, on a stalled plan.
  - No within-cap case was infeasible.
- **Change.** In `evaluateConstraints`, for a segment whose start knot pins a velocity, that axis's cap becomes max(SegMaxVel, |v_BC|, |v_plan(start of segment)|) (`Guidance/PolynomialQP/evaluateConstraints.m:37-41`, `:75-87`; `evaluateConstraints.cpp:58-62`, `:96-108`). The last term is the D3 amendment, recorded under D3: with the FIX-F03 BC, the plan's CG velocity at t = 0 is v + ρ(ḃ3 − ḃ3,plan(0)), not exactly v_BC. Without the amendment, 22 ASTRA and 21 TOAD tilted starts were infeasible from t = 0, all velocity-dominated and pinned. Acceleration is not changed: it is free at knot 1, so dilation can reduce it.
- **Acceptance** (`fix_regressions` F01):
  - Frozen: ASTRA 82 → 0 infeasible, TOAD 93 → 0.
  - New within-cap, over-cap and tilted sets (100 each per vehicle): 0 infeasible.
  - Zero remaining cases attributable to the over-cap artifact. MATLAB/C++ status and `Feasible` identical on all 1,000 cases.
  - No other infeasible case remains, so none needed a physical-limit classification. Outside the sampled ranges, the constructed start z0 = 0.6 m, vz = −8 m/s is still infeasible, as it should be.
  - `verify_polyqp` now reports feasibility separately for within-cap and over-cap starts.

### FIX-F03: knot-1 P-space BC from the actual attitude

- **Cause.** `solveAllLegs` set P = r_cg + ρ e3 at every knot and applied the velocity BC to dP/dt, assuming b3 = e3 and ḃ3 = 0 at the abort instant.
- **Change.**
  - Optional `B3BC` and `B3DotBC` on `TOADWaypoint`, default NaN = current behaviour (`TOADWaypoint.m:14-15`, `:57-58`; `PolyQP_types.h:166-167`, `:180-181`).
  - Carried per knot in `Plan.B3BC`/`B3DotBC` (`resolveWaypoints.m:55-71`, `:166-167`, `:185-186`, `:206-207`, `:258-259`; `resolveWaypoints.cpp:67-89`, `:197-198`, `:220-221`, `:244-245`, `:296-297`).
  - `solveAllLegs` sets P = r_cg + ρ b3 and dP/dt = v_cg + ρ ḃ3 where given, and the refinement loop leaves a given attitude alone (`PolyQP.m:143-155`, `:224`; `PolyQP.cpp:169-185`, `:270-272`).
  - The abort block computes b3 = C_IB e3 with C_IB = `quatRot(q)'`, and ḃ3 = (C_IB ω_B) × b3 (`TrajectoryManager.m:31-36`; `AbortGuidance.cpp:31-41`). The convention was checked, not guessed: `quatRot` returns C_BI (`FlightEstimator.m:20` uses `R_b2i = quatRot(q)'`), and `DCM_Quat_Conversion(quatRot(q)') = q`. The C++ uses a `Real` transpile of `Helper/quatRot.m` (`quatRot.cpp`, flag F36). A non-finite quaternion or rate degrades to the old BC.
- **Acceptance.** As agreed after Phase 1, the prompt's criterion (plan CG at t = 0 = `X_est` to ≤ 1e-9) is achievable only with the knot-1 acceleration and jerk pinned. Results (`fix_regressions` F03):
  - **(i)** The P-space BC is exact (0 error, MATLAB and C++) on 200 tilted, rotating starts.
  - **(ii)** With the FIX-F02 pin on, plan CG position and velocity at t = 0 equal `X_est` to ≤ 3.6e-15, tilted and upright, MATLAB and C++. Without the pin, see D6.
  - **(iii)** Upright, zero-rate starts are bit-identical to the pre-F03 BC in 400 / 400 frozen cases.
  - MATLAB/C++ ≤ 6.3e-10.

### FIX-F02 (DECISION D4): optional acceleration/jerk pin, off

- **Cause (Phase 1.4).** Knot-1 acceleration is free, so the plan's t = 0 thrust and attitude differ from the vehicle's. Upright frozen: attitude step median 8.2° / max 25.5° (ASTRA) and 2.9° / 17.3° (TOAD); thrust step up to 41 % / 17 % of hover, which takes 0.16 s / 0.26 s to slew at the planner's own thrust-rate limit.
- **Change.** `ABORT_PIN_THRUST` (MATLAB constant `false`, C++ parameter defaulting to `false`, flag F37). When on, waypoint 1 gets `AccelBC` = (T_hover/m) b3 − g e3 = g(b3 − e3) and `JerkBC` = g ḃ3 (`TrajectoryManager.m:16`, `:37-40`; `AbortGuidance.cpp:15`, `:43-47`). The jerk pin is the addition agreed after Phase 1; the velocity match needs it.
- **Acceptance** (`fix_regressions` F02). With the pin on, the attitude step is 0°, the thrust-direction step ≤ 2.2° and the thrust step ≤ 5.4 % of hover. MATLAB/C++ parity through `AbortGuidance` is ≤ 9e-9.

### FIX-F29: TOAD plans from the estimated propellant

- **Cause.** `buildVehicleLimits` used full tanks for `m_init`, ρ and the mass integration's start, so a mid-flight TOAD abort planned with full-tank mass and inertia.
- **Change.** `Opts.InitOxMass`/`InitFuMass` (default: full tank) set `m_lox0`/`m_ipa0`, and therefore `m_init`, ρ and `flatnessInvert`'s `ml_s(1)`, `mi_s(1)`. `Constants.OxMass`/`FuMass` stay the full-tank reference for `ComputeJtot` (unchanged). The abort block passes `X_est(14)`, `X_est(15)`.
  - `buildVehicleLimits.m:1`, `:21-27`; `PolyQP.m:16`; `TrajectoryManager.m:26-27`
  - `buildVehicleLimits.cpp:7-39`; `PolyQP.cpp:25`; `AbortGuidance.cpp:22-26`; `PolyQP_types.h:150-151` (flag F38)
- **Acceptance** (`fix_regressions` F29):
  - TOAD at 50 %, 20 %, 5 % and 0 % tanks: `States(1, 14:15) = X_est(14:15)` exactly in MATLAB and C++, mass non-increasing, `m_init` and ρ from the estimate, MATLAB/C++ ≤ 3.2e-9.
  - Full-tank TOAD and ASTRA (which ignores the estimate) are bit-identical in 400 / 400.
  - Finding: depleted TOAD aborts are much longer, by a median 29–37 s and up to 406 s. Near-empty TOAD hovers at 1,778 N, barely above the planner's minimum thrust T_lo = 0.6 × 2,936 = 1,761 N, so it can hardly descend. All were feasible.

### FIX-F27: consistent heading fallback in `flatnessInvert`

- **Cause.** When b3 ∥ x_C, b2 switches to `cross(u, x_C_alt)`, but its derivatives still used `x_C`.
- **Change.** `x_ref` is the heading actually used, in `cd` and `cdd` (`flatnessInvert.m:50-66`; `flatnessInvert.cpp:71-88`).
- **Acceptance** (`fix_regressions` F27). On a synthetic leg whose specific force is exactly along x_C at t0, ω and ω̇ match finite differences of the fallback-frame attitude to 2.0e-11 and 6.4e-7 (pre-fix: 83 % error). "Continuous across the switch" cannot hold, because the switch changes the heading reference and the attitude itself jumps by up to 90° about b3. The test is that inside the branch the rates are the derivatives of the attitude the branch defines. All 400 frozen plans are bit-identical to the pre-fix function.

## Decisions (defaults chosen; the lead may reverse them)

| ID | Default chosen | Evidence | How to reverse |
| --- | --- | --- | --- |
| D1 | Fallback when PolyQP fails: hold the current lateral position and descend rest-to-rest from z0 to 0 with the r = 5 (degree-9) profile, all derivatives zero at both ends. T = max(0.2, c_v·\|z0\|/1 m/s), c_v = 315/128, using the waypoint-3 vertical cap. Hover thrust m·g, zero gimbal and roll, upright, propellant held at the estimate. Non-finite inputs are replaced (position 0, propellant full tank) so the output stays finite. | F08 acceptance. The profile starts at rest, so the vehicle's current velocity is a step the controller must absorb. 40 m → 98 s. | Edit `abortFallback.m` and mirror it in `abortFallback.cpp`. |
| D2 | `Feasible = 0` plans are flown (`INFEASIBLE_FLOWN`), not replaced by the fallback. | After FIX-F01 no sampled start is infeasible; the constructed z0 = 0.6 m, vz = −8 m/s is. | Change `if ~planOk` to `if ~planOk \|\| ~Report.Feasible` in `TrajectoryManager.m:67` and mirror it at `AbortGuidance.cpp:91`. |
| D3 | A pinned start-velocity segment's cap = max(SegMaxVel, \|v_BC\|, \|v_plan(segment start)\|) per axis. The last term is the amendment for the F03 BC. Not applied to acceleration. | F01 acceptance; the 22/21 tilted cases above. | Restore `vCap = Plan.SegMaxVel(:, m)` in `evaluateConstraints.m:76` and `.cpp:97`; drop the `vStart` term to keep only the prompt's rule. |
| D4 | `ABORT_PIN_THRUST = false`. This reverses my Phase 1 recommendation: with the pin on, 87 / 100 ASTRA and 99 / 100 TOAD tilted starts are infeasible (velocity utilisation up to 12 and 48). A pinned lateral start acceleration with a free landing x, y makes the minimum-crackle plan drift sideways, and stretching time makes it worse: 19 m of drift at the nominal times, 304 m at 4× for a 15° ASTRA start. Upright starts are unaffected (0 infeasible) because the pin is then a = 0, j = 0. | F02 / F03 tables. | Set `ABORT_PIN_THRUST = true` (`TrajectoryManager.m:16`) and pass `true` to the C++ `AbortGuidance`. A planner change (for example pinning the landing x, y, or excluding pinned-acceleration legs from dilation) would be needed first. |
| D5 | `ABORT_HOVER_ALT = 2` m (no value change; the lead's spec says 1 m). | F01 flag. | Change the constant in `TrajectoryManager.m:15` and `AbortGuidance.cpp:14`. |
| D6 | F03 is passed to waypoint 1 as the prompt specifies, although the pin is off. The plan's P (centre of percussion) then starts at the vehicle's actual P. The plan's CG at t = 0 differs from `X_est` by ρ\|b3_actual − b3_plan(0)\|, which on tilted starts is **larger** than before F03: max 0.16 m / 0.20 m/s (ASTRA), 0.34 m / 0.56 m/s (TOAD), against 0.09 / 0.06 and 0.15 / 0.13 before. Upright starts are unchanged. | F03 table. | Remove `'B3BC', b3_est, 'B3DotBC', b3dot_est` from waypoint 1 (`TrajectoryManager.m:36`, and leave `options1.B3BC`/`B3DotBC` at NaN in `AbortGuidance.cpp`). That restores the old BC exactly. |
| D7 | In Simulink, the MATLAB block falls back only for a non-finite estimate; a PolyQP error stops the simulation (no `try/catch` in code generation). The flight C++ falls back on every error. | Measured `coder.screener`/`codegen` on R2026b; `check_abort_codegen.m`. | — (a code-generation limit; converting PolyQP's `error()` calls into status returns would be a larger MATLAB change). |
| D8 | FIX-F07 assembles the reduced system segment by segment instead of Tzᵀ R_PP Tz. | F07 section: 1.4e-1 error against `vpa` with the literal formula. | Not recommended. |

## Blockers

### B1: roll-channel units (F30), controller side

The quantity is a roll **torque in N·m** at every computational stage:
- planner: `allocateControls.m:24`, `:36`, from the Euler moment at `flatnessInvert.m:130`;
- abort hand-off: `TrajectoryManager.m`, passed through unchanged;
- LQR design model: `Controls/TOAD_Controller_Gen.m:61`, `:67`, where roll enters `MB`;
- controller arithmetic: `cpp/ASTRAv2_Controller.cpp:96`, `:109`, `:185-191`, where `U_ff(3)` is summed with `U_dist(2) = ang_accel_dist(2) · J(2)`, a torque;
- actuator: `Actuator_Models_2.slx`, `Actuators`, a unit-free lag;
- plant: `Flight Dynamics/FlightDynamicsGen.m:10` ("roll input (Nm)"), `:107`, `:109`, `:114`.

J_zz = 0.020 kg·m² (ASTRA) and 20.6 kg·m² (TOAD, full tanks; 20.0 empty). Planned |U(4)| in the abort cases: max 0.008 N·m (ASTRA), 0.23 N·m (TOAD).

Wrong (controller side, not edited):
1. `cpp/ASTRAv2_Controller.cpp:191` labels the channel `// Roll (rad/s^2)` and clamps it at ±4. That is ±4 N·m: 200 rad/s² on ASTRA and 0.19 rad/s² on TOAD. If ±4 rad/s² was intended, the clamp must be ±4·J_zz (0.08 N·m ASTRA, ≈ 82 N·m TOAD); otherwise only the label is wrong.
2. `Controls/TOAD_Controller_Gen.m:8` comment `% roll input (rad/sec^2)` → N·m.
3. `cpp/README.md:67` ("roll accel") and `:99` (`roll_rad_sec_squared`: "Roll acceleration command, rad/s^2") → torque, N·m.

Related, inside PolyQP but a design value, not a unit bug: `buildVehicleLimits.m:8`, `:46` builds `RollLimit = 0.9 * MaxRollRate` (4.5 N·m) and `RollRateLim = 5` N·m/s from a constant named as a rate. 4.5 exceeds the controller's ±4 clamp. Its limits should come from the actuator and J_zz once item 1 is settled. Not changed.

## Other findings (not changed)

- **F06:** `StraightLine` on abort waypoint 3 is still inert (NaN lateral positions). FIX-F07 removes its consequence for conditioning; with all derivatives zero at both ends of segment 2, the optimum is Δ = 0, so the vehicle lands below the ABT GND point.
- **F04, F05:** no change (redundant passes 2–3; no time contraction).
- **Memory:** see the table above and `TRANSPILE_FLAGS.md`.
- **Long depleted-tank TOAD aborts** (FIX-F29) also hit the F11 grid cap (3999 samples), so the evaluation grid gets coarser above 200 s.

## Not verified

- Cortex-M7 compile: `arm-none-eabi-g++` is not installed, so step 3 of `build_host.sh` was skipped.
- No closed-loop Simulink simulation of an abort was run. The fixes are verified at the planner level (MATLAB vs C++ and the acceptance criteria). How the controller responds to the remaining t = 0 steps, the fallback's velocity step or the D6 CG offset was not measured.
- `check_abort_codegen.m` runs MATLAB Coder `codegen` on `TrajectoryManager` (lib, code only) with a synthetic `constantsTOAD.Traj`. It does not build `TOAD_TrajectoryGen_3.slx` itself.
- No MC_v5 dispersion data is in the repo, so the tilt and rate distributions are sampled (tilt ≤ 30°, |ω| ≤ 1 rad/s).
- The flight trigger (`ABORT = 0`, F09/F31) is unchanged and still never fires.

## Reproducing

From `TOAD/cpp/PolyQP/verify`:

```
matlab -batch "build_polyqp_mex; build_polyqp_mex('float'); verify_polyqp; check_abort_codegen; diagnose_aborts('fixed')"
bash build_host.sh
```

- `verify_polyqp` runs every transpile check on the fixed code plus all `fix_regressions` sections. It writes `verify_polyqp_report.md` and errors on any failure.
- `diagnose_aborts('fixed')` adds the "after" rows to `diagnosis_report.md`; the baseline rows come from `results/diag_baseline.mat`.
- The frozen baseline is `results/baseline_report.md`, `results/baseline_run.log` and `results/baseline_cases_slim.mat`. The full `baseline_cases.mat` (33 MB) is not committed; to regenerate it, run `verify_polyqp('baseline')` on commit `a46d227`.
