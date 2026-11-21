// Abort: flight harness for the PolyQP abort guidance (ASTRA). This is the only header the flight software includes.
// See INTEGRATION.md. Not reentrant: plan() uses one static solver workspace (PolyQP/TRANSPILE_FLAGS.md).
#pragma once

#include "matlab_funcs.h"
#include <stdint.h>

namespace Abort {

// Output node capacity: rows the caller's trajectory arrays must hold for a planned abort (0.1 s grid). Sized from
// ABORT_MAX_START_ALT_M (PolyQP/PolyQP_config.h); Abort.cpp asserts it equals the solver's POLYQP_MAX_NODES.
constexpr int MAX_POINTS = 601;

enum class Status : uint8_t {
  OK = 0,               // planned and feasible
  INFEASIBLE_FLOWN = 1, // planned, but a limit is exceeded somewhere along the plan; the plan is written anyway (D2)
  FALLBACK = 2          // the solver could not plan (or capacity < MAX_POINTS): closed-form vertical descent (D1)
};

struct Result {
  Status status;
  int num_points;        // rows written into the caller's arrays
  float duration_s;      // time[num_points - 1]
  uint8_t solver_status; // PolyQPStatus (0 = Ok; 8 = CapacityExceeded), for logging
};

// Plans an abort from the current state and writes it into the caller's trajectory arrays (time starts at 0). Never fails:
// if the solver cannot produce a plan, the arrays hold the closed-form vertical-descent fallback. The arrays have the layout
// of ASTRA_Trajectory_Data.h (TRAJ_TIME[N], TRAJ_STATES[N][15], TRAJ_INPUTS[N][4]) and must be writable RAM with
// capacity >= MAX_POINTS rows; a smaller capacity gets the fallback, sampled to fit (at least 2 rows are needed; with fewer,
// nothing is written and num_points = 0). Blocks until done (milliseconds on a host; measure on the target).
Result plan(const Vector15 &x_est, // [q(4); pos(3); vel(3); omega(3); m_lox; m_ipa], as fed to ASTRAv2_Controller
            const Vector4 &u_cmd,  // last command [gimbal pitch (rad); gimbal yaw (rad); thrust (N); roll], for the accel BC
            const Vector3 &a_dist, // estimated disturbance accel (LESO_Position), or zero
            float time[], float states[][15], float inputs[][4], int capacity);

} // namespace Abort
