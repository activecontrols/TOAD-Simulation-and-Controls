// Transpiled from TOAD/Guidance/AbortGuidance.m (the abort block formerly in TrajectoryManager.m) @ 5ac6bd23, with the abort fixes of ABORT_FIXES.md
#include "AbortGuidance.h"
#include "PolyQP.h"

namespace PolyQP {

AbortStatus AbortGuidance(const Vector15d &X_est, const Vector4d &u_cmd, const Vector3d &a_dist, const NodeSink &Sink, ReportT &Report, PolyQPStatus &planStatus) {
    // Abort guidance: plans the abort from the estimated state X_est and returns the trajectory (Time, States, Inputs) and its
    // status (0 = OK, 1 = INFEASIBLE_FLOWN, 2 = FALLBACK). The caller decides when to abort and calls this once.
    // This is the abort block of Guidance/TrajectoryManager.m from branch abort-fixes (TrajectoryManager.m was later replaced by
    // the interpolating trajectory manager, commit fd08005), with the start-of-abort acceleration / jerk boundary condition
    // from attitude and thrust (FIX-F02, ASTRA integration Phase 5). The MATLAB-simulation trigger machinery (ABORT, tAbort,
    // abortCalculated and its persistence) was removed in Phase 6 because the caller decides when to abort.
    // It is the MATLAB reference that cpp/PolyQP/AbortGuidance.cpp transpiles.
    // u_cmd is the last command [gimbal pitch; gimbal yaw; thrust (N); roll] and a_dist the estimated disturbance acceleration
    // (LESO_Position); they set the start acceleration and jerk (default: hover thrust with the gimbal centred, no disturbance).
    // Report is PolyQP's report when a plan was made (empty for the fallback).
    // TRANSPILE-FLAG F09 / F31 (removed, Phase 6): the trigger (ABORT, tAbort, abortCalculated and the FIX-F32 one-attempt latch) is gone; the caller decides when to abort (Abort::plan).
    constexpr double ABORT_HOVER_ALT = 2; // DECISION-D5 (TRANSPILE-FLAG F01): abort waypoint 2 altitude (m); the lead's spec says 1 m above the pad
    static_assert(ABORT_HOVER_ALT == POLYQP_ABORT_HOVER_ALT_M, "STATIC: the grid capacity derivation (PolyQP_config.h) assumes this altitude");
    // FIX-F02: default command: hover thrust, gimbal centred, no roll
    // FIX-F02: no disturbance estimate
    AbortStatus abortStatus = AbortStatus::OK; // FIX-F08: 0 = OK, 1 = INFEASIBLE_FLOWN (D2), 2 = FALLBACK (D1)

    // ASTRA-ONLY: the MATLAB Opts (R = 5, InitOxMass/InitFuMass from X_est(14:15), FIX-F29) is not carried: r is hardcoded and ASTRA carries no propellant.
    const TOADWaypointT dummyWP = TOADWaypoint(Vector3d(0, 0, 0));
    TOADWaypointArray &Waypoints = workspace().AbortGuidance.Waypoints; // STATIC: shared workspace
    bool planOk = repmat(dummyWP, 1, 3, Waypoints);

    // FIX-F02: body thrust vector from the last command, as LESO_Position computes it: T * [cos(th) sin(ph); -sin(th); cos(th) cos(ph)]
    const Vector4d q_est = X_est(R1(1, 4)) / mnorm(Vector4d(X_est(R1(1, 4))));
    const Matrix3_3d C_IB = quatRot(q_est).transpose();
    const Vector3d T_B = V1(u_cmd, 3) * Vector3d(std::cos(V1(u_cmd, 1)) * std::sin(V1(u_cmd, 2)), -std::sin(V1(u_cmd, 1)), std::cos(V1(u_cmd, 1)) * std::cos(V1(u_cmd, 2)));
    const Vector3d w_I = C_IB * Vector3d(X_est(R1(11, 13)));
    // FIX-F02: start acceleration a0 = C_IB * T_B / m + g + a_dist (g = [0; 0; -g], z up; m = m_dry, ASTRA carries no propellant)
    const Vector3d C_IB_T_B = C_IB * T_B;
    const Vector3d a0 = C_IB_T_B / ASTRA::m_dry + Vector3d(0, 0, -ASTRA::g) + a_dist;
    // FIX-F02: start jerk j0 = (w_I x C_IB * T_B) / m with w_I = C_IB * w_B: the thrust is held constant over the abort instant (DECISION-D6)
    const Vector3d j0 = mcross(w_I, C_IB_T_B) / ASTRA::m_dry;
    // FIX-F03: thrust axis b3 = C_IB * T_B / |T_B| (the body z axis when the gimbal is centred) and db3/dt = w_I x b3, consistent with a0 and j0
    Vector3d b3_est;
    if (mnorm(T_B) > 0) {
        b3_est = C_IB_T_B / mnorm(T_B);
    } else {
        b3_est = C_IB * Vector3d(0, 0, 1);
    }
    const Vector3d b3dot_est = mcross(w_I, b3_est);
    // DECISION-D9: the heading reference is the estimated heading, x_C = C_IB(:, 2) x e3, so the plan starts at the estimated attitude
    Vector3d hdg_est = mcross(Vector3d(Col1(C_IB, 2)), Vector3d(0, 0, 1));
    if (mnorm(hdg_est) < double(1e-6)) {
        hdg_est = Vector3d(1, 0, 0);
    }
    TOADWaypointOptions options1;
    options1.VelBC = X_est(R1(8, 10));
    options1.VelFixed.setConstant(true);
    options1.AccelBC = a0;
    options1.AccelFixed.setConstant(true);
    options1.JerkBC = j0;
    options1.JerkFixed.setConstant(true);
    options1.Label = "ABT START";
    options1.B3BC = b3_est;
    options1.B3Fixed.setConstant(true);
    options1.B3DotBC = b3dot_est;
    options1.B3DotFixed.setConstant(true);
    options1.HDGRef = hdg_est;
    Waypoints(1) = TOADWaypoint(X_est(R1(5, 7)), options1); // FIX-F02/F03: knot-1 state, acceleration, jerk and attitude
    // TRANSPILE-FLAG F02 (fixed by FIX-F02, ASTRA integration Phase 5): the knot-1 acceleration and jerk were free (or optionally pinned to hover, ABORT_PIN_THRUST, now removed); they are the vehicle's current acceleration and jerk.
    // FIX-F02: the end knots keep AccelBC = 0 and JerkBC = 0: the a0/j0 formula at hover (C_IB = I, T_B = m g e3, a_dist = 0, w = 0) gives exactly these
    // TRANSPILE-FLAG F01: abort waypoint 2 altitude is 2 m (ABORT_HOVER_ALT, DECISION-D5); the lead's spec says 1 m above the pad. Transpiled as-is per MATLAB.
    TOADWaypointOptions options2;
    options2.Label = "ABT GND";
    options2.HoldTime = double(0.0);
    options2.HoldTimeGiven = true;
    options2.VelBC = Vector3d(0, 0, 0);
    options2.VelFixed.setConstant(true);
    options2.AccelBC = Vector3d(0, 0, 0);
    options2.AccelFixed.setConstant(true);
    options2.JerkBC = Vector3d(0, 0, 0);
    options2.JerkFixed.setConstant(true);
    options2.HDGRef = hdg_est;
    Waypoints(2) = TOADWaypoint(Vector3d(0, 0, ABORT_HOVER_ALT), options2, Vector3b(false, false, true)); // NO-NAN: [nan; nan; ABORT_HOVER_ALT], the lateral landing position is free
    TOADWaypointOptions options3;
    options3.Label = "ABT END";
    options3.HoldTime = double(0.0);
    options3.HoldTimeGiven = true;
    options3.MaxVel = Vector3d(3, 3, 1);
    options3.StraightLine = true;
    options3.VelBC = Vector3d(0, 0, 0);
    options3.VelFixed.setConstant(true);
    options3.AccelBC = Vector3d(0, 0, 0);
    options3.AccelFixed.setConstant(true);
    options3.JerkBC = Vector3d(0, 0, 0);
    options3.JerkFixed.setConstant(true);
    options3.HDGRef = hdg_est;
    Waypoints(3) = TOADWaypoint(Vector3d(0, 0, 0), options3, Vector3b(false, false, true)); // NO-NAN: [nan; nan; 0], the lateral landing position is free
    // NO-NAN: the MATLAB "non-finite start state -> fallback" branch (FIX-F08) was removed: the estimate is always real
    // TRANSPILE-FLAG F08 (fixed by FIX-F08): MATLAB ignored Report.Feasible and PolyQP errors and had no fallback.
    // TRANSPILE-FLAG F16: the PolyQPStatus (planStatus) is returned to the C++ caller for logging. Transpiled as-is per MATLAB.
    planStatus = PolyQPStatus::Ok;
    if (planOk) {
        InitialGuessT InitialGuess;
        QPTrajT &QPTraj = workspace().AbortGuidance.QPTraj; // STATIC: shared workspace
        // FIX-F08: code generation has no try/catch, so a PolyQP error is caught only in interpreted MATLAB
        // TRANSPILE-FLAG F35: the C++ catches every PolyQP error (including CapacityExceeded) through the status return, like interpreted MATLAB's try/catch.
        // STREAM: Time = InitialGuess.Time(:), States = InitialGuess.X' and Inputs = [U'; U(:, end)'] go straight to Sink row by row
        // STREAM: (sampleOutputs); a failed plan writes nothing, so the fallback below writes the whole trajectory.
        planStatus = PolyQP(Waypoints, Sink, InitialGuess, QPTraj, Report);
        planOk = planStatus == PolyQPStatus::Ok;
        if (planOk) {
            // FIX-F08: an infeasible plan is flown anyway (D2) and flagged
            if (Report.Feasible) {
                abortStatus = AbortStatus::OK;
            } else {
                abortStatus = AbortStatus::INFEASIBLE_FLOWN;
            }
        }
    }
    if (!planOk) { // FIX-F08: closed-form vertical descent that cannot fail (D1)
        abortFallback(X_est, V1(Waypoints(3).MaxVel, 3), Sink);
        abortStatus = AbortStatus::FALLBACK;
    }
    return abortStatus;
}

}
