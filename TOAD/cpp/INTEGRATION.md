# Abort integration (ASTRA flight software)

> **BLOCKER B1: roll-channel units.** `inputs[i][3]` is PolyQP's `U(4)`, a roll **torque in N·m**. The follower feeds it to `ci.U_ff[3]`, and `ASTRAv2_Controller.cpp:191` labels that channel rad/s² and clamps it at ±4. The nominal trajectories carry the same quantity: `ASTRA_Circle_v003` is also a PolyQP output, with column 4 = roll. So the abort is consistent with what the follower already flies, and `Abort::plan` outputs what PolyQP produces. Which unit the controller should take is still unresolved (`PolyQP/ABORT_FIXES.md`, B1).
>
> If the controller wants roll acceleration, the conversion is α_z = τ_z / J_zz = τ_z / 0.02. Apply it once, where the follower sets `ci.U_ff[3]`, so the nominal and abort trajectories are converted alike. Also fix the ±4 clamp. Do **not** convert inside `Abort.cpp` alone: the abort would then disagree with the nominal trajectories.

## The call

```cpp
Vector4 u_prev(last_co.gimbal_yaw_deg * DEG_TO_RAD, last_co.gimbal_pitch_deg * DEG_TO_RAD, last_co.thrust_N, last_co.roll_rad_sec_squared);
Abort::Result r = Abort::plan(ControllerAndEstimator::get_controller_x_est(), u_prev, ControllerAndEstimator::get_position_disturbance(), ABORT_TIME, ABORT_STATES, ABORT_INPUTS, Abort::MAX_POINTS);
traj_time = ABORT_TIME; traj_states = ABORT_STATES; traj_inputs = ABORT_INPUTS; traj_num_points = r.num_points; traj_t0 = last_time_s; pt = 0;
```

`examples/TrajectoryFollower_abort_example.cpp` shows the whole follower; every changed line is marked `ABORT:`.

- **`plan()` never fails.**
  - `r.status`: 0 = `OK`, 1 = `INFEASIBLE_FLOWN` (planned but a limit is exceeded; it is flown anyway, D2), 2 = `FALLBACK` (closed-form vertical descent, D1).
  - `r.solver_status`: the `PolyQPStatus` (0 = Ok, 8 = CapacityExceeded). Log both.
- **`u_cmd`** is the previous tick's command in the order `LESO_Position` uses: `[U(0); U(1); thrust N; roll]`.
  - In this repo, `controller_and_estimator.cpp` writes `U(0)` into `co.gimbal_yaw_deg` and `U(1)` into `co.gimbal_pitch_deg`, so the two fields are swapped relative to `ASTRAv2_Controller.cpp:188-189`, which calls `U(0)` pitch. Hence the order above. **Check this against the FC copy.**
  - `co` is also EMA-filtered (α = 0.5 per tick), so `u_prev` lags `LESO`'s `last_U_cmd` by up to one tick.
  - Thrust is assumed constant over the abort instant (D6).
- **`a_dist`** is the LESO_Position disturbance estimate. A zero vector is also valid.

## Getting `x_est` and `a_dist`

`controller_and_estimator.*` gained two read-only getters. They copy values the existing code already computes and change nothing else:

- `const Vector15 &get_controller_x_est()` returns `X_ctrl` exactly as passed to `ASTRAv2_Controller` on the last call: `[q(4); pos(3); vel(3); bias-corrected body rates(3); 0; 0]`.
- `const Vector3 &get_position_disturbance()` returns the `LESO_Position` `a_dist` of the last call. It is the same value as `co.dist_accel_x/y/z`.

Both are zero before the first `get_controller_output` call.

**The FC's `controller_and_estimator.*` is a different version from this repo's:** the follower there passes `ci.X_trg` / `ci.U_ff`, while here they are arguments. Port the two getters by hand: two globals, two assignments after `Dist_LESO << a_dist, u_dist;`, and two one-line functions.

**Frames.**
- **Position and velocity:** north / west / up (`est_pos_north/west/up`, matching `ci.gps.pos`), z up. This is the trajectory frame: the controller compares `X_ctrl` with `X_trg` directly, and the abort uses g = [0; 0; −g].
- **Attitude:** scalar-first quaternion with `quatRot(q) = C_BI`, as in `LESO_Position.cpp`.
- **Body rates:** gyro after `DigitalNF`, minus the estimated gyro bias `x_est(10:12)`, in the order `[gyro_yaw, gyro_pitch, gyro_roll]` = body x, y, z. Body z is the thrust / roll axis, J_zz = 0.02.

## Files to copy into the FC build

| File | Note |
|---|---|
| `Abort.h`, `Abort.cpp` | The only header the follower includes. |
| `matlab_funcs.h`, `matlab_helpers.cpp` | Replace the FC copies. Additions only: the existing float helpers are unchanged, and the double overloads of `quatRot` / `DCM_Quat_Conversion` only take double arguments. |
| `PolyQP/*.cpp` (14 files) | `PolyQP, resolveWaypoints, assembleQP, solveLeg, allocateTimes, evaluateConstraints, flatnessInvert, allocateControls, sampleOutputs, PolyEval, buildVehicleLimits, TOADWaypoint, AbortGuidance, abortFallback` |
| `PolyQP/*.h` | `PolyQP.h`, `PolyQP_types.h`, `PolyQP_config.h`, `AbortGuidance.h`, and the constants header `ASTRA_PolyQP_Constants.h`. |
| `controller_and_estimator.*` | Port the two getters only (above). |

- Include paths: the directory holding `matlab_funcs.h`, which includes `PolyQP/PolyQP_config.h` for the capacity typedefs, and Eigen 3.4 through `ArduinoEigenDense.h`, as the flight code already does.
- Nothing under `PolyQP/verify/` or `verify/` goes into the FC build.

**Compiler flags.**
- Build the PolyQP sources and `Abort.cpp` with `-ffp-contract=off`. GCC defaults to `fast` for GNU C++, and on the M7's double-precision FPU that fuses multiply-adds, which changes rounding away from the verified results.
- **No `-ffast-math`** (no `-ffinite-math-only`). The compile-time guard was removed with the NaN checks, so the build itself must ensure this. Check the Teensyduino *Optimize* menu.
- Exceptions and RTTI are not used, and no heap is used: the host test runs under `-DEIGEN_NO_MALLOC` and counts 0 allocations inside `plan()`.

**RAM.**
- Solver static `.bss` + `.data`: **25,056 B**. The solver workspace in `PolyQP.o` is 24,880 B, `Abort.o` 144 B and `assembleQP.o` 32 B. Measured with `size` on the host build (`PolyQP/verify/build_host.sh`, step 2b); the target is ≤ 256 KB.
- Abort trajectory buffer, option (b) below: `Abort::MAX_POINTS` × 80 B = **48,080 B**.
- Stack: about **13 KB** at the deepest call chain, measured with host clang `-fstack-usage` (plan → AbortGuidance → PolyQP → solveAllLegs → allocateTimes → solveLeg → Eigen). Measure it on the target.
- `plan()` is not reentrant; it uses one static workspace.

## Writable trajectory arrays

`plan()` writes into arrays with the `ASTRA_Trajectory_Data.h` layout: `time[N]`, `states[N][15]`, `inputs[N][4]`, with N ≥ `Abort::MAX_POINTS` = 601 for a planned abort. With a smaller `capacity` it writes the fallback, sampled to fit.

- **(a) Overwrite in place.** The FC's `Trajectory::TRAJ_*` arrays must be non-`const` RAM arrays with ≥ 601 rows, and the nominal trajectory is copied into them at boot. This costs the same 48 KB, and the nominal trajectory is lost after an abort.
- **(b) Recommended.** Keep the nominal arrays `const` in flash and declare one static abort buffer (`ABORT_TIME[Abort::MAX_POINTS]`, `ABORT_STATES[...][15]`, `ABORT_INPUTS[...][4]`). The follower reads through pointers (`traj_time`, `traj_states`, `traj_inputs`, `traj_num_points`) and switches them on abort. The example does this.

The abort does not produce the gain-schedule arrays of `TrajectoryDataset` (`k_trans`, `k_rot`, `la`, `lt`). If the FC follower schedules gains from the trajectory, it must hold the last nominal gains, or constants, during an abort.

## Time base and segment index

The abort's `time[0] = 0` is the tick whose estimate it planned from, so set `traj_t0 = last_time_s`. Look up the trajectory with `t - traj_t0`, and reset `pt = 0`.

The original nested `for (pt) / while` cannot do this, so the example uses one loop. At the top of each tick, `pt` advances while `last_time_s - traj_t0 >= traj_time[pt + 1]`, and the flight ends when `pt` reaches the last point of the trajectory being flown.

## Solve time: `plan()` blocks the 1 kHz loop

**Host (Apple M5 Pro, `-O2`, flight flags):** median **1.0 ms**, p95 5.4 ms, worst **13.7 ms** over the 240 verification starts (200 upright, 40 tilted; `PolyQP/verify/build_host.sh` step 4).

**The Teensy worst case must be measured.** A 600 MHz Cortex-M7 is expected to be one to two orders of magnitude slower, so the worst case may reach hundreds of milliseconds. Long plans dilate up to 40 outer iterations.

While `plan()` runs:
- no estimator or controller update happens;
- the props and gimbal hold their last command;
- the next tick sees a large `loop_dT`.

After the solve, the follower resumes partway into the plan, at the elapsed time. The example also resets `lastloop` so that the loop does not run un-throttled to catch up. If the measured worst case is too long, the options are a lower `ABORT_MAX_START_ALT_M`, a cap on the outer iterations, or running `plan()` off the control loop. Each of those is a design decision.

## Decisions for the lead

- **g mismatch:** `ASTRA_PolyQP_Constants.h` uses g = **9.80145** (`LoadTOADParams.m`). The controller (`init_controller_and_estimator_constants`), `TrajectoryManager.cpp` via the trajectory header metadata, and `verify_trajectory_manager.m` use **9.8015**. m_dry = 1.1690, rTB = 0.26, J = diag(0.067, 0.067, 0.02) and MaxThrust = 16 agree everywhere. Not resolved here.
- **`ABORT_MAX_START_ALT_M = 40` → `MAX_POINTS = 601`:**
  - Worst-case T_total = 2.5 × T_init(40 m) = 59 s (`PolyQP/PolyQP_config.h`).
  - The 2.5 factor is measured: with the Phase 5 boundary condition, the worst observed dilation was 2.25 × T_init, and the worst plan was 48.3 s.
  - A longer plan returns CapacityExceeded, and the abort falls back.
  - The lead sets the altitude; the grid, the static_asserts and `MAX_POINTS` follow from it.
- **D6:** thrust is held constant over the abort instant (Ṫ = 0) for the start jerk.
- **D9:** the abort's heading reference is the estimated heading, so the plan starts at the estimated attitude.
- **B1:** roll units, above.

The planner itself produces some start offsets; none of them is an error in `Abort`:

- The plan's yaw rate at t = 0 is not the estimate, because the heading is held constant. Roll and pitch rates match.
- With a deflected gimbal, the plan's start attitude aligns body z with the thrust axis, so it differs from `x_est` by the gimbal angle.
- With a nonzero `a_dist`, the plan's start thrust axis is tilted by about m|a_dist|/|T|, and its CG by about 0.22 m times that angle. This is because the planner has no disturbance model.
- Most tilted starts are `INFEASIBLE_FLOWN`: with the start acceleration pinned and the landing x, y free, the plan drifts sideways (Phase 5 finding). Fixing that is a planner change.
- The landing x, y is free, so the vehicle lands wherever the minimum-crackle plan ends. That can be tens of metres from the start for a fast lateral start.

## Verification

- `PolyQP/verify/build_host.sh` runs:
  - the host smoke test and the harness test `verify/abort_harness_test.cpp`, both with `-Werror -DEIGEN_NO_MALLOC`;
  - the flight-flags `-Werror` build and the RAM table;
  - the Cortex-M7 compile, if `arm-none-eabi-g++` is installed;
  - the host timing.
- MATLAB parity of the solver and the abort block: `PolyQP/verify/verify_polyqp.m`, report `verify_polyqp_report.md`.
- The example follower is not compiled in this repo; it depends on FC headers.
