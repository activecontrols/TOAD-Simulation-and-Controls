// Transpiled from TOAD/Guidance/AbortGuidance.m (the abort block formerly in TrajectoryManager.m) @ 5ac6bd23, with the abort fixes of ABORT_FIXES.md
#pragma once

#include "PolyQP_types.h"

namespace PolyQP {

// FIX-F08: outcome of the (single) abort plan; the MATLAB persistent abortStatus (0, 1, 2)
enum class AbortStatus : int {
    OK = 0,               // PolyQP succeeded and Report.Feasible
    INFEASIBLE_FLOWN = 1, // PolyQP succeeded, Report.Feasible = 0; the plan is flown anyway (D2)
    FALLBACK = 2          // PolyQP failed (planStatus, e.g. CapacityExceeded); closed-form vertical descent (D1)
};

struct AbortGuidanceState {
    bool abortCalculated = false;
    double tAbort = 0;
    AbortStatus abortStatus = AbortStatus::OK;  // FIX-F08: valid once abortCalculated
    PolyQPStatus planStatus = PolyQPStatus::Ok; // FIX-F08: PolyQP's status on the abort attempt, for logging (C++ only)
};

// FIX-F08: the trajectory (Time, States, Inputs) is written to Sink exactly once, when the abort fires (plan or fallback); the
// return value is state.abortStatus. STREAM: Sink receives row k = (Time(k), States(k, :), Inputs(k, :)); STATIC: no array is
// stored here. FIX-F02: u_cmd (last command [gimbal pitch; gimbal yaw; thrust N; roll]) and a_dist (estimated disturbance
// acceleration) set the start acceleration and jerk; the defaults are hover thrust with the gimbal centred and no disturbance.
AbortStatus AbortGuidance(double t, const Vector15d &X_est, AbortGuidanceState &state, const NodeSink &Sink, ReportT &Report, double ABORT = 0,
                          const Vector4d &u_cmd = Vector4d(0, 0, ASTRA::m_dry * ASTRA::g, 0), const Vector3d &a_dist = Vector3d::Zero());

// FIX-F08: the closed-form fallback, written to Sink (at most Sink.capacity rows: the grid coarsens to fit, F34).
void abortFallback(const Vector15d &X_est, double vCapZ, const NodeSink &Sink);

}
