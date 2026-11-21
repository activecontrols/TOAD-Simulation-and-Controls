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

// FIX-F08: plans the abort once from X_est (the caller decides when; the MATLAB trigger machinery is gone) and writes the
// trajectory (Time, States, Inputs) to Sink: STREAM: Sink receives row k = (Time(k), States(k, :), Inputs(k, :)); STATIC: no
// array is stored here. FIX-F02: u_cmd (last command [gimbal pitch; gimbal yaw; thrust N; roll]) and a_dist (estimated
// disturbance acceleration) set the start acceleration and jerk. Report is valid when a plan was made; planStatus is PolyQP's
// status (Ok, or why the fallback was flown, e.g. CapacityExceeded).
AbortStatus AbortGuidance(const Vector15d &X_est, const Vector4d &u_cmd, const Vector3d &a_dist, const NodeSink &Sink, ReportT &Report, PolyQPStatus &planStatus);

// FIX-F08: the closed-form fallback, written to Sink (at most Sink.capacity rows: the grid coarsens to fit, F34).
void abortFallback(const Vector15d &X_est, double vCapZ, const NodeSink &Sink);

}
