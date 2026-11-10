# ASTRA Controller and Estimator: MATLAB MEX Wrapper

This directory provides a MATLAB interface to the flight C++ estimator and controller. It allows you to evaluate the flight implementation alongside simulations without uploading code to the vehicle.

## What is MEX?

A MEX function lets MATLAB call compiled C++ as if it were a regular MATLAB function. Compile the source to create a binary (on Windows, `astra_mex.mexw64`); MATLAB can then execute the C++ implementation through calls to `astra_mex(...)`.

The wrapper uses the vehicle's C++ implementation, including its single-precision `float` calculations, so tests exercise the same software logic used in flight.

## What is in this folder

| File | What it is |
|---|---|
| `astra_mex.cpp` | The MATLAB entry point; translates MATLAB structs and vectors into C++ types and converts returned values back. This is the only MATLAB-specific C++ source file. |
| `controller_and_estimator.cpp / .h` | Main flight pipeline, executing the notch filter, estimator, disturbance observers, and controller on each iteration. |
| `GroundEstimator.cpp` | Kalman filter used during ground operation, combining IMU, magnetometer, and GPS data. |
| `FlightEstimator.cpp` | Kalman filter used during flight, using IMU and GPS measurements. |
| `DigitalNF.cpp` | Digital notch and low-pass filtering of IMU measurements. |
| `Q_gen.cpp` | Generates the estimator's process-noise matrix. |
| `ASTRAv2_Controller.cpp` | Flight controller with an outer position loop and an inner attitude loop. |
| `LESO_Position.cpp`, `LESO_Attitude.cpp` | Position and attitude disturbance observers. |
| `matlab_helpers.cpp` | Utility functions for quaternion and rotation calculations. |
| `matlab_funcs.h` | Shared declarations and Eigen type definitions used by the source files. |
| `astra_structs.h` | Data structures for inputs, outputs, and estimator state. |
| `Extras/ArduinoEigenDense.h` | A compatibility header that maps the Arduino-style `#include <ArduinoEigenDense.h>` include to standard Eigen on a PC. It also declares `Matrix9_4`. |
| `Extras/eigen-3.4.0/` | Eigen 3.4.0 linear algebra headers; no separate installation step is needed. |
| `build_mex.m` | Build script that compiles the source files into the MEX binary. |
| `test_astra.m` | Basic test that supplies synthetic sensor readings to the MEX function and plots its outputs. |

## Prerequisites

- MATLAB (Windows instructions below; macOS/Linux work with a gcc/clang compiler)
- A MATLAB-supported C++ compiler. Check with `mex -setup C++`. (This project was built with Microsoft Visual C++.)
- Eigen 3.4.0, placed at `Extras/eigen-3.4.0/`. The folder that directly contains `Eigen/` must be the one the build script points at. If the unzip created a doubled `eigen-3.4.0/eigen-3.4.0/`, the build script uses the inner one.

## Step 1: Compile the MEX wrapper

1. In MATLAB, change the working directory to the folder containing `astra_mex.cpp`:
   ```matlab
   cd 'C:\path\to\this\folder'
   ```
2. Execute the build script:
   ```matlab
   build_mex
   ```
3. A successful build prints `MEX completed successfully.` and creates `astra_mex.mexw64` in this directory. Warnings about converting `double` values to `float` are expected; focus on messages marked as errors.

If MATLAB reports that the MEX binary is locked during a rebuild, run `clear mex` and try again.

## Step 2: Call the MEX function

The MEX has two commands: `'init'` and `'step'`.

### `astra_mex('init')`

Initializes the estimator, filters, and controller. Run this once at the start of each simulation.

### `[co, internals] = astra_mex('step', ci, X_trg, U_ff, ideal_dT, loop_dT)`

Executes one control-loop iteration. Call it at every simulation time step. The MEX retains estimates, filter history, and integrator state between calls.

**Inputs**

- `ci`: struct of sensor data and status flags (see below)
- `X_trg`: 15x1 double target state, ordered `[quaternion(4); position(3); velocity(3); body rates(3); 2 unused]`. The quaternion is `[w x y z]`; use `[1 0 0 0]` for upright. An all-zero quaternion is invalid once in flight.
- `U_ff`: 4x1 double feedforward `[gimbal pitch; gimbal yaw; thrust N; roll accel]`
- `ideal_dT`, `loop_dT`: time step in seconds. Use `0.001` (1 kHz). The noise matrices and notch filter assume this rate.

**`ci` fields** (all numeric scalars; use 1/0 or true/false for flags)

```matlab
ci.GND_val          % true = on the ground (ground estimator runs, controller outputs 0)
ci.new_mag_packet   % true on steps where a new magnetometer reading arrived
ci.new_gps_packet   % true on steps where a new GPS reading arrived

ci.imu.accel_x, accel_y, accel_z          % m/s^2
ci.imu.gyro_yaw, gyro_pitch, gyro_roll    % rad/s
ci.mag.mag_x, mag_y, mag_z

ci.gps.pos.north, west, up                % m
ci.gps.vel.north, west, up                % m/s
ci.gps.posCovNN, posCovNE, posCovND, posCovEE, posCovED, posCovDD
ci.gps.velCovNN, velCovNE, velCovND, velCovEE, velCovED, velCovDD

% optional (default 0, not used by the current code):
ci.target_pos_north, target_pos_west, target_pos_up
```

If a required field is absent, MATLAB reports the name of the missing field.

**Outputs**

`co` (controller output):

| Field | Meaning |
|---|---|
| `thrust_N` | Thrust command, N |
| `roll_rad_sec_squared` | Roll acceleration command, rad/s^2 |
| `gimbal_pitch_deg` | Gimbal pitch command, deg (low-pass filtered) |
| `gimbal_yaw_deg` | Gimbal yaw command, deg (low-pass filtered) |

`internals` (optional second output, for debugging/plots):

- `internals.x_est`: struct with the estimated state: quaternion (`q_vec_w/x/y/z`), position (`est_pos_north/west/up`), velocity (`est_vel_north/west/up`), gyro, accel and mag biases
- `internals.filter_out`: 9x1 filtered IMU data
- `internals.VelTarget`, `AccelTarget`, `TargetAtt`, `VelErrorI`, `AttErrorI`: controller internals (these are filled in only if the controller code writes them)

### Minimal working example

```matlab
clear astra_mex            % start from a clean state
astra_mex('init');

ci.GND_val = true; ci.new_mag_packet = false; ci.new_gps_packet = false;
ci.imu = struct('accel_x',0,'accel_y',0,'accel_z',9.8, ...
                'gyro_yaw',0,'gyro_pitch',0,'gyro_roll',0);
ci.mag = struct('mag_x',0.385,'mag_y',0.032,'mag_z',-0.922);
ci.gps.pos = struct('north',0,'west',0,'up',0);
ci.gps.vel = struct('north',0,'west',0,'up',0);
for f = ["posCovNN","posCovNE","posCovND","posCovEE","posCovED","posCovDD", ...
         "velCovNN","velCovNE","velCovND","velCovEE","velCovED","velCovDD"]
    ci.gps.(f) = 1;
end

X_trg = zeros(15,1); X_trg(1) = 1;
U_ff  = [0; 0; 11.46; 0];

[co, internals] = astra_mex('step', ci, X_trg, U_ff, 0.001, 0.001);
```

### Calling it from a simulation loop

Initialize with `astra_mex('init')` before entering the loop. At each time step, populate `ci` with the simulated sensor values and call `astra_mex('step', ...)`. Apply the returned `co` commands to the vehicle dynamics model.

## Step 3: Execute the basic test

```matlab
test_astra
```

The test simulates five seconds at 1 kHz with the vehicle held upright and stationary: two seconds in ground mode followed by three seconds in flight mode. It prints basic checks and creates four plots showing the estimated quaternion, estimated position, thrust, and gimbal commands. A dashed line indicates the mode transition. This verifies basic execution and output behavior; it does not validate a flight trajectory.

## Important behavior and notes

- **Ground and flight modes:** When `ci.GND_val = true`, the ground estimator runs and controller outputs are set to zero. When it is `false`, the flight estimator and controller run. The transition copies the estimator covariance and resets the controller integrators.
- **Resetting a simulation:** Before a new run, execute `clear astra_mex` (or `clear mex`) and then call `astra_mex('init')` to clear the state retained by the C++ implementation.
- **GPS offset in flight mode:** The flight estimator compensates for a 0.31 m GPS lever arm along the body z-axis. For an upright vehicle whose origin is at zero, the GPS up-position should therefore be approximately 0.31 m.
- **Numeric precision:** Calculations in C++ use `float`, while returned values are represented as MATLAB `double`. Small differences from a fully double-precision MATLAB simulation are expected.
- **Recompile after source changes:** MATLAB does not rebuild the binary automatically. Run `build_mex` again after editing C++ files.

## Compatibility edits for the MATLAB build

These adjustments accommodate Microsoft's compiler (MSVC), which is used by the Windows MATLAB build instead of GCC. They do not affect the Teensy build.

1. In `astra_structs.h`, define `__packed` differently for MSVC because the GCC attribute is not supported there:
   ```cpp
   #ifdef _MSC_VER
   #define __packed
   #else
   #define __packed __attribute__((__packed__))
   #endif
   ```
2. In `ASTRAv2_Controller.cpp`, an `Axis` expression used a ternary operator with incompatible Eigen expression types. MSVC rejects that combination, so it was replaced with an equivalent `if` statement:
   ```cpp
   Vector3 Axis = Vector3::Zero();
   if (n > 1e-6f) {
     Axis = att_vec / n;
   }
   ```
3. `Extras/ArduinoEigenDense.h` provides the compatibility layer that includes standard Eigen and declares `Matrix9_4`.

Keep these three compatibility changes when modifying the Teensy-side code; removing them can break the MEX build.

## Troubleshooting guide

| Problem | Fix |
|---|---|
| `astra_mex.cpp not detected` | You're in the wrong folder. `cd` into this folder and check with `pwd` and `dir *.cpp`. |
| `Cannot open include file: 'Eigen/Dense'` | The Eigen path in `build_mex.m` is wrong. Make sure `Extras/eigen-3.4.0/eigen-3.4.0/Eigen/Dense` exists (`isfile(...)` returns 1), or adjust the path. |
| Errors about `__attribute__` or `__packed__` | The `astra_structs.h` change above is missing. |
| `C2446` / "no conversion from ... CwiseNullaryOp" | A `? :` mixing Eigen expressions. Rewrite it as an `if`/`else`. |
| `Missing struct field '...'` at runtime | `ci` is missing a field. Add it (see the field list above). |
| NaNs in the output | Check `X_trg(1)` is 1 (valid quaternion), `dT` is `0.001`, and that sensor inputs are finite. |
| Can't rebuild: file in use | Run `clear mex`, then `build_mex`. |
