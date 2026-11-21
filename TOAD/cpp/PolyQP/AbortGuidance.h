// Transpiled from TOAD/Guidance/AbortGuidance.m (the abort block formerly in TrajectoryManager.m) @ 5ac6bd23, with the abort fixes of ABORT_FIXES.md
#pragma once

#include "PolyQP_types.h"
#include "matlab_funcs.h"

namespace PolyQP {

using AbortXEstT = VecX<15>;
using AbortTimeT = VecX<POLYQP_MAX_NODES>;
using AbortStatesT = MatX<POLYQP_MAX_NODES, 15>;
using AbortInputsT = MatX<POLYQP_MAX_NODES, 4>;

using AbortTimeF = Eigen::Matrix<float, Eigen::Dynamic, 1, Eigen::ColMajor, POLYQP_MAX_NODES, 1>;
using AbortStatesF = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, POLYQP_MAX_NODES, 15>;
using AbortInputsF = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, POLYQP_MAX_NODES, 4>;

// FIX-F08: outcome of the (single) abort plan; the MATLAB persistent abortStatus (0, 1, 2)
enum class AbortStatus : int {
    OK = 0,               // PolyQP succeeded and Report.Feasible
    INFEASIBLE_FLOWN = 1, // PolyQP succeeded, Report.Feasible = 0; the plan is flown anyway (D2)
    FALLBACK = 2          // PolyQP failed (planStatus) or the start state is non-finite; closed-form vertical descent (D1)
};

struct AbortGuidanceState {
    bool abortCalculated = false;
    Real tAbort = 0;
    AbortStatus abortStatus = AbortStatus::OK;     // FIX-F08: valid once abortCalculated
    PolyQPStatus planStatus = PolyQPStatus::Ok;     // FIX-F08: PolyQP's status on the abort attempt, for logging (C++ only)
};

// FIX-F08: Time/States/Inputs are written exactly once, when the abort fires (plan or fallback); the return value is state.abortStatus.
// FIX-F02: ABORT_PIN_THRUST (default false, the MATLAB constant) pins the knot-1 accel/jerk to hover thrust along the current b3 (DECISION-D4).
AbortStatus AbortGuidance(Real t, const AbortXEstT &X_est, const PolyQPConstants &constantsTOAD, AbortGuidanceState &state, AbortTimeT &Time,
                          AbortStatesT &States, AbortInputsT &Inputs, ReportT &Report, Real ABORT = 0, bool ABORT_PIN_THRUST = false);

AbortStatus AbortGuidance(float t, const Vector15 &X_est, const PolyQPConstants &constantsTOAD, AbortGuidanceState &state, AbortTimeF &Time,
                          AbortStatesF &States, AbortInputsF &Inputs, ReportT &Report, float ABORT = 0, bool ABORT_PIN_THRUST = false);

void abortFallback(const AbortXEstT &X_est, const PolyQPConstants &Constants, Real vCapZ, AbortTimeT &Time, AbortStatesT &States, AbortInputsT &Inputs);

}
