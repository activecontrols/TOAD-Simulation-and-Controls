// Abort: flight harness for the PolyQP abort guidance (ASTRA). See Abort.h and INTEGRATION.md.
// The solver runs in double; float appears only here, at the boundary (the caller's flight types and trajectory arrays).
#include "Abort.h"
#include "PolyQP/AbortGuidance.h"
#include "PolyQP/PolyQP.h"

namespace Abort {

static_assert(MAX_POINTS == POLYQP_MAX_NODES, "Abort::MAX_POINTS must equal the solver's output node capacity (PolyQP_config.h)");
static_assert(static_cast<int>(Status::OK) == static_cast<int>(PolyQP::AbortStatus::OK) &&
                  static_cast<int>(Status::INFEASIBLE_FLOWN) == static_cast<int>(PolyQP::AbortStatus::INFEASIBLE_FLOWN) &&
                  static_cast<int>(Status::FALLBACK) == static_cast<int>(PolyQP::AbortStatus::FALLBACK),
              "Abort::Status mirrors PolyQP::AbortStatus");
static_assert(static_cast<int>(PolyQP::PolyQPStatus::CapacityExceeded) == 8, "Abort.h documents solver_status 8 = CapacityExceeded");

namespace {

// Vertical speed cap of the fallback descent: ABT END MaxVel(3) in the abort block (Guidance/AbortGuidance.m, waypoint 3),
// which PolyQP::AbortGuidance passes to abortFallback. Used here only when the caller's capacity is below MAX_POINTS.
constexpr double FALLBACK_VCAP_Z = 1.0;

// Row writer handed to the solver (PolyQP::NodeSink): output node k (1-based) goes straight into row k - 1 of the caller's
// arrays, converted to float, so no intermediate copy of the trajectory exists.
struct FloatRows {
  float *time;
  float (*states)[15];
  float (*inputs)[4];
  int rows;
};

void put_row(void *ctx, int k, double t, const Vector15d &X, const Vector4d &U) {
  FloatRows *out = static_cast<FloatRows *>(ctx);
  const int i = k - 1;
  out->time[i] = static_cast<float>(t);
  for (int j = 0; j < 15; ++j) {
    out->states[i][j] = static_cast<float>(X(j));
  }
  for (int j = 0; j < 4; ++j) {
    out->inputs[i][j] = static_cast<float>(U(j));
  }
  if (k > out->rows) {
    out->rows = k;
  }
}

} // namespace

Result plan(const Vector15 &x_est, const Vector4 &u_cmd, const Vector3 &a_dist, float time[], float states[][15], float inputs[][4], int capacity) {
  const Vector15d x_est_d = x_est.cast<double>();
  const Vector4d u_cmd_d = u_cmd.cast<double>();
  const Vector3d a_dist_d = a_dist.cast<double>();

  FloatRows out{time, states, inputs, 0};
  Result result{Status::FALLBACK, 0, 0.0f, static_cast<uint8_t>(PolyQP::PolyQPStatus::CapacityExceeded)};

  if (capacity < 2) {
    return result; // nothing fits: the fallback needs a start and an end row
  }

  if (capacity < MAX_POINTS) {
    // A planned abort may need up to MAX_POINTS rows, so with less room the fallback is written, sampled to fit (F34).
    const PolyQP::NodeSink sink{&out, capacity, &put_row};
    PolyQP::abortFallback(x_est_d, FALLBACK_VCAP_Z, sink);
  } else {
    const PolyQP::NodeSink sink{&out, MAX_POINTS, &put_row};
    static PolyQP::ReportT report; // STATIC: kept off the stack; not reentrant (one abort plan at a time)
    PolyQP::PolyQPStatus solver_status = PolyQP::PolyQPStatus::Ok;
    const PolyQP::AbortStatus status = PolyQP::AbortGuidance(x_est_d, u_cmd_d, a_dist_d, sink, report, solver_status);
    result.status = static_cast<Status>(static_cast<uint8_t>(status));
    result.solver_status = static_cast<uint8_t>(solver_status);
  }

  result.num_points = out.rows;
  result.duration_s = out.rows > 0 ? time[out.rows - 1] : 0.0f;
  return result;
}

} // namespace Abort
