// Verification only (never part of the flight build): host test of the flight harness Abort::plan (Abort.h).
// Built and run by PolyQP/verify/build_host.sh with -Werror -DEIGEN_NO_MALLOC (any Eigen heap allocation asserts), and
// global operator new / delete replaced to count every heap allocation made inside Abort::plan (must be 0).
//
// Checks, for representative states written into arrays shaped like ASTRA_Trajectory_Data.h:
//   time starts at 0 and is strictly increasing; num_points <= capacity; duration_s = time[num_points - 1];
//   the first row matches x_est (planned paths: attitude, position, velocity, roll/pitch rate; fallback: position);
//   the last row is at rest at z = 0; all three status paths occur (the fallback both from an undersized capacity and from
//   the solver's CapacityExceeded); and a mocked follower call sequence (abort trigger, pointer switch, re-based time,
//   time-driven segment index) flies the abort buffer to its end.
#include "Abort.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>

// ------------------------------------------------------------------ heap allocation counter
static long g_allocations = 0;
void *operator new(std::size_t n) {
  ++g_allocations;
  void *p = std::malloc(n == 0 ? 1 : n);
  if (p == nullptr) {
    std::abort();
  }
  return p;
}
void *operator new[](std::size_t n) { return operator new(n); }
void operator delete(void *p) noexcept { std::free(p); }
void operator delete[](void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }
void operator delete[](void *p, std::size_t) noexcept { std::free(p); }

// ------------------------------------------------------------------ arrays shaped like ASTRA_Trajectory_Data.h
// Option (b) of INTEGRATION.md: one static abort buffer in RAM next to the const nominal arrays.
static float ABORT_TIME[Abort::MAX_POINTS];
static float ABORT_STATES[Abort::MAX_POINTS][15];
static float ABORT_INPUTS[Abort::MAX_POINTS][4];

// A const nominal trajectory (stands in for TRAJ_ASTRA_CIRCLE_V003_*): hover at 5 m, 71 points over 7 s.
constexpr int NOMINAL_NUM_POINTS = 71;
static float NOMINAL_TIME[NOMINAL_NUM_POINTS];
static float NOMINAL_STATES[NOMINAL_NUM_POINTS][15];
static float NOMINAL_INPUTS[NOMINAL_NUM_POINTS][4];

static int g_failures = 0;
static void check(bool ok, const char *what, const char *tag) {
  if (!ok) {
    std::printf("  FAIL [%s] %s\n", tag, what);
    ++g_failures;
  }
}

static const char *status_name(Abort::Status s) {
  switch (s) {
  case Abort::Status::OK: return "OK";
  case Abort::Status::INFEASIBLE_FLOWN: return "INFEASIBLE_FLOWN";
  case Abort::Status::FALLBACK: return "FALLBACK";
  }
  return "?";
}

static Vector15 make_state(float roll_deg, float pitch_deg, float yaw_deg, float x, float y, float z, float vx, float vy, float vz, float wx, float wy,
                           float wz) {
  // q = qz(yaw) * qy(pitch) * qx(roll), scalar first (the estimator's x_est(0:3) convention)
  const float d2r = 3.14159265358979f / 180.0f;
  const float cr = std::cos(0.5f * roll_deg * d2r), sr = std::sin(0.5f * roll_deg * d2r);
  const float cp = std::cos(0.5f * pitch_deg * d2r), sp = std::sin(0.5f * pitch_deg * d2r);
  const float cy = std::cos(0.5f * yaw_deg * d2r), sy = std::sin(0.5f * yaw_deg * d2r);
  Vector15 X;
  X << cr * cp * cy + sr * sp * sy, sr * cp * cy - cr * sp * sy, cr * sp * cy + sr * cp * sy, cr * cp * sy - sr * sp * cy, x, y, z, vx, vy, vz, wx, wy, wz, 0, 0;
  return X;
}

struct Case {
  const char *tag;
  Vector15 x_est;
  Vector4 u_cmd;
  Vector3 a_dist;
  int capacity;
  Abort::Status expect;
  int expect_solver; // -1: any
};

// Checks the written rows; returns the result for the summary.
static Abort::Result run_case(const Case &c) {
  const long alloc0 = g_allocations;
  const Abort::Result r = Abort::plan(c.x_est, c.u_cmd, c.a_dist, ABORT_TIME, ABORT_STATES, ABORT_INPUTS, c.capacity);
  const long allocs = g_allocations - alloc0;
  std::printf("%-44s status=%-16s solver=%u points=%3d/%3d duration=%7.3f s heap=%ld\n", c.tag, status_name(r.status), r.solver_status, r.num_points, c.capacity,
              r.duration_s, allocs);
  const int n = r.num_points;

  check(r.status == c.expect, "status as expected", c.tag);
  check(c.expect_solver < 0 || r.solver_status == c.expect_solver, "solver_status as expected", c.tag);
  check(allocs == 0, "no heap allocation inside Abort::plan", c.tag);
  check(n >= 2 && n <= c.capacity && n <= Abort::MAX_POINTS, "2 <= num_points <= capacity", c.tag);
  if (n < 2) {
    return r;
  }
  check(ABORT_TIME[0] == 0.0f, "time starts at 0", c.tag);
  bool increasing = true;
  for (int i = 1; i < n; ++i) {
    increasing = increasing && ABORT_TIME[i] > ABORT_TIME[i - 1];
  }
  check(increasing, "time strictly increasing", c.tag);
  check(r.duration_s == ABORT_TIME[n - 1], "duration_s = time[num_points - 1]", c.tag);

  // First row = x_est. Position always (the fallback holds the estimated position); attitude, velocity and the roll/pitch
  // rates for a planned path (the fallback starts upright at rest: the D1 velocity step). The yaw rate is not matched: the
  // abort holds a constant heading (Phase 5 finding).
  float e_pos = 0, e_vel = 0, e_q = 0, e_w = 0;
  for (int j = 0; j < 3; ++j) {
    e_pos = std::fmax(e_pos, std::fabs(ABORT_STATES[0][4 + j] - c.x_est(4 + j)));
    e_vel = std::fmax(e_vel, std::fabs(ABORT_STATES[0][7 + j] - c.x_est(7 + j)));
  }
  for (int j = 0; j < 4; ++j) {
    e_q = std::fmax(e_q, std::fabs(ABORT_STATES[0][j] - c.x_est(j)));
  }
  for (int j = 0; j < 2; ++j) {
    e_w = std::fmax(e_w, std::fabs(ABORT_STATES[0][10 + j] - c.x_est(10 + j)));
  }
  // A disturbance estimate tilts the plan's start thrust axis away from the vehicle's (the planner has no disturbance model:
  // a0 includes a_dist, so the plan's thrust is C_IB T_B + m a_dist), by at most alpha = m |a_dist| / |T_B|. The plan's start
  // attitude and its CG (P-space offset rho b3, rho ~ 0.22 m) then differ from x_est by O(alpha): same exclusion as set S4 of
  // the BC regression. Likewise a deflected gimbal: the flat-output model puts the thrust on the body z axis, so the plan's
  // start attitude aligns body z with the thrust axis b3 = C_IB T_B / |T_B|, a rotation by the gimbal angle
  // gamma = acos(cos(th) cos(ph)) away from x_est (set S3). With neither, the match is exact to float rounding.
  const float alpha = 1.1690f * c.a_dist.norm() / std::fmax(c.u_cmd(2), 1e-3f);
  const float gamma = std::acos(std::fmin(1.0f, std::cos(c.u_cmd(0)) * std::cos(c.u_cmd(1))));
  const float tol_q = 1e-6f + 0.5f * (alpha + gamma) * 1.01f, tol_p = 0.25f * alpha;
  check(e_pos <= 1e-5f * (1.0f + std::fabs(c.x_est(6))) + (r.status == Abort::Status::FALLBACK ? 0.0f : tol_p), "first row position = x_est", c.tag);
  if (r.status != Abort::Status::FALLBACK) {
    check(e_q <= tol_q, "first row attitude = x_est", c.tag);
    check(e_vel <= 1e-5f + tol_p, "first row velocity = x_est", c.tag);
    const bool centred = c.u_cmd(0) == 0.0f && c.u_cmd(1) == 0.0f;
    check(!centred || e_w <= 1e-5f + tol_p, "first row roll/pitch rate = x_est (gimbal centred)", c.tag);
  }

  // Last row at rest at z = 0
  float v_end = 0, w_end = 0;
  for (int j = 0; j < 3; ++j) {
    v_end = std::fmax(v_end, std::fabs(ABORT_STATES[n - 1][7 + j]));
    w_end = std::fmax(w_end, std::fabs(ABORT_STATES[n - 1][10 + j]));
  }
  check(std::fabs(ABORT_STATES[n - 1][6]) <= 1e-4f, "last row z = 0", c.tag);
  check(v_end <= 1e-4f, "last row velocity = 0", c.tag);
  check(w_end <= 1e-3f, "last row body rate = 0", c.tag);
  std::printf("    first-row err: pos %.2g vel %.2g q %.2g w(roll/pitch) %.2g | last row: z %.2g |v| %.2g |w| %.2g\n", e_pos, e_vel, e_q, e_w,
              ABORT_STATES[n - 1][6], v_end, w_end);
  return r;
}

// ------------------------------------------------------------------ mocked follower call sequence
// Stands in for ControllerAndEstimator (get_controller_x_est / get_position_disturbance) and the follower loop of
// examples/TrajectoryFollower_abort_example.cpp: a single loop that advances pt by time, so pt and the time base can be
// reset when the abort fires.
namespace MockController {
Vector15 x_est;
Vector3 a_dist;
const Vector15 &get_controller_x_est() { return x_est; }
const Vector3 &get_position_disturbance() { return a_dist; }
} // namespace MockController

static float lerp(float a, float b, float t0, float t1, float t) {
  if (t <= t0) return a;
  if (t >= t1) return b;
  if (t0 == t1) return b;
  return a + (b - a) * ((t - t0) / (t1 - t0));
}

static void mocked_follower() {
  const char *tag = "mocked follower";
  for (int i = 0; i < NOMINAL_NUM_POINTS; ++i) {
    NOMINAL_TIME[i] = 0.1f * static_cast<float>(i);
    for (int j = 0; j < 15; ++j) NOMINAL_STATES[i][j] = 0.0f;
    NOMINAL_STATES[i][0] = 1.0f;
    NOMINAL_STATES[i][6] = 5.0f;
    for (int j = 0; j < 4; ++j) NOMINAL_INPUTS[i][j] = 0.0f;
    NOMINAL_INPUTS[i][2] = 11.46f;
  }
  // Read pointers (option b): start on the const nominal arrays
  const float *traj_time = NOMINAL_TIME;
  const float (*traj_states)[15] = NOMINAL_STATES;
  const float (*traj_inputs)[4] = NOMINAL_INPUTS;
  int traj_num_points = NOMINAL_NUM_POINTS;

  const float dt = 0.001f; // 1 kHz loop
  const float abort_at = 3.2f;
  float t0 = 0.0f; // time base of the trajectory being flown
  int pt = 0;
  bool auto_land = false;
  Abort::Result r{};
  float last_U_ff[4] = {0, 0, 0, 0};
  float last_X_trg[15] = {0};
  Vector4 last_co(0.0f, 0.0f, 11.46f, 0.0f);
  long ticks = 0;
  float traj_t_end = 0.0f;
  for (float now = 0.0f;; now += dt, ++ticks) {
    if (!auto_land && now >= abort_at) {
      // abort trigger: plan from the estimate and the previous tick's command, switch the read pointers, re-base time
      MockController::x_est = make_state(4, -3, 20, 1.5f, -2.0f, 5.0f, 0.4f, -0.3f, -0.2f, 0.02f, -0.01f, 0.0f);
      MockController::a_dist = Vector3(0.05f, -0.02f, 0.1f);
      r = Abort::plan(MockController::get_controller_x_est(), last_co, MockController::get_position_disturbance(), ABORT_TIME, ABORT_STATES, ABORT_INPUTS,
                      Abort::MAX_POINTS);
      traj_time = ABORT_TIME;
      traj_states = ABORT_STATES;
      traj_inputs = ABORT_INPUTS;
      traj_num_points = r.num_points;
      t0 = now;
      pt = 0;
      auto_land = true;
      std::printf("%-44s status=%-16s solver=%u points=%3d duration=%7.3f s (fired at t = %.3f s)\n", tag, status_name(r.status), r.solver_status, r.num_points,
                  r.duration_s, now);
    }
    const float t = now - t0;
    while (pt < traj_num_points - 2 && t >= traj_time[pt + 1]) {
      ++pt; // advance the segment index by time
    }
    if (t >= traj_time[traj_num_points - 1]) {
      traj_t_end = t;
      break; // end of the trajectory being flown
    }
    for (int i = 0; i < 15; ++i) {
      last_X_trg[i] = lerp(traj_states[pt][i], traj_states[pt + 1][i], traj_time[pt], traj_time[pt + 1], t);
    }
    for (int i = 0; i < 4; ++i) {
      last_U_ff[i] = lerp(traj_inputs[pt][i], traj_inputs[pt + 1][i], traj_time[pt], traj_time[pt + 1], t);
    }
    if (ticks > 200000) break;
  }
  check(auto_land, "abort fired", tag);
  check(traj_num_points == r.num_points && pt == r.num_points - 2, "segment index walked to the last abort segment", tag);
  check(std::fabs(traj_t_end - r.duration_s) <= 2 * dt, "flew the abort for its full duration on the re-based clock", tag);
  check(std::fabs(last_U_ff[2] - traj_inputs[r.num_points - 1][2]) <= 0.05f, "final feedforward thrust = last abort row", tag);
  check(std::fabs(last_X_trg[6]) <= 1e-3f, "final target altitude = 0 (landed)", tag);
  std::printf("    flew %ld ticks; pt = %d of %d; abort clock at end %.3f s\n", ticks, pt, r.num_points - 2, traj_t_end);
}

int main() {
  const float m_g = 1.1690f * 9.80145f;
  const Vector4 hover(0.0f, 0.0f, m_g, 0.0f);
  const Vector3 zero3 = Vector3::Zero();
  const Case cases[] = {
      {"upright, 10 m, hover command", make_state(0, 0, 0, 3, -2, 10, 1, 0.5f, -1, 0, 0, 0), hover, zero3, Abort::MAX_POINTS, Abort::Status::OK, 0},
      {"upright, 30 m, moving, yawed 40 deg", make_state(0, 0, 40, -15, 12, 30, -3, 4, 2, 0, 0, 0), hover, zero3, Abort::MAX_POINTS, Abort::Status::OK, 0},
      {"slight tilt, rates, disturbance", make_state(2, -1.5f, -60, 0, 0, 8, 0.3f, -0.2f, 0.5f, 0.05f, -0.03f, 0.1f), hover, Vector3(0.05f, -0.02f, 0.1f),
       Abort::MAX_POINTS, Abort::Status::OK, 0},
      {"upright, 0.6 m, vz = -8 m/s", make_state(0, 0, 0, 2, -3, 0.6f, 0, 0, -8, 0, 0, 0), hover, zero3, Abort::MAX_POINTS, Abort::Status::OK, 0},
      {"tilted 12/-8 deg, rotating (infeasible, D2)", make_state(12, -8, 30, 4, -6, 15, 1, -1, 0, 0.3f, -0.2f, 0.1f), hover, zero3, Abort::MAX_POINTS,
       Abort::Status::INFEASIBLE_FLOWN, 0},
      {"tilted, gimbal 3/-2 deg, 0.9 m g thrust", make_state(-6, 10, -120, -8, 3, 25, -2, 0.5f, -1, -0.1f, 0.2f, 0), Vector4(0.0524f, -0.0349f, 0.9f * m_g, 0.0f),
       zero3, Abort::MAX_POINTS, Abort::Status::INFEASIBLE_FLOWN, -1},
      {"undersized capacity (71 rows)", make_state(0, 0, 0, 3, -2, 10, 1, 0.5f, -1, 0, 0, 0), hover, zero3, NOMINAL_NUM_POINTS, Abort::Status::FALLBACK, 8},
      {"300 m start (solver CapacityExceeded)", make_state(0, 0, 0, 1, 2, 300, 0, 0, 0, 0, 0, 0), hover, zero3, Abort::MAX_POINTS, Abort::Status::FALLBACK, 8},
  };
  int seen[3] = {0, 0, 0};
  for (const Case &c : cases) {
    const Abort::Result r = run_case(c);
    ++seen[static_cast<int>(r.status)];
  }
  check(seen[0] > 0 && seen[1] > 0 && seen[2] > 0, "all three status paths exercised", "summary");

  // Capacity below 2 rows: nothing written, reported as fallback with num_points = 0
  const Abort::Result r1 = Abort::plan(cases[0].x_est, hover, zero3, ABORT_TIME, ABORT_STATES, ABORT_INPUTS, 1);
  check(r1.status == Abort::Status::FALLBACK && r1.num_points == 0, "capacity 1: FALLBACK with num_points = 0", "capacity 1");

  mocked_follower();

#if defined(EIGEN_NO_MALLOC)
  std::printf("built with EIGEN_NO_MALLOC\n");
#endif
  if (g_failures == 0) {
    std::printf("abort_harness_test: PASS\n");
    return 0;
  }
  std::printf("abort_harness_test: %d FAILURES\n", g_failures);
  return 1;
}
