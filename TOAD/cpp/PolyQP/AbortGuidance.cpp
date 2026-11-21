// Transpiled from TOAD/Guidance/AbortGuidance.m (the abort block formerly in TrajectoryManager.m) @ 5ac6bd23, with the abort fixes of ABORT_FIXES.md
#include "AbortGuidance.h"
#include "PolyQP.h"

namespace PolyQP {

AbortStatus AbortGuidance(double t, const Vector15d &X_est, AbortGuidanceState &state, const NodeSink &Sink, ReportT &Report, double ABORT, bool ABORT_PIN_THRUST) {
    // TRANSPILE-FLAG F09: MATLAB hard-codes the trigger ABORT = 0 (never fires). Transpiled as-is per MATLAB.
    // TRANSPILE-FLAG F31: ABORT is a parameter defaulting to 0 (the MATLAB value) so flight code cannot abort unless a caller sets it; tests set it to 1.
    // TRANSPILE-FLAG F31 (cont.): the persistent tAbort/abortCalculated/abortStatus live in AbortGuidanceState (tAbort = 0, abortCalculated = false, abortStatus = OK). Transpiled as-is per MATLAB.
    // TRANSPILE-FLAG F37: ABORT_PIN_THRUST is a parameter defaulting to the MATLAB constant (false), like ABORT (F31), so tests can switch it on.
    constexpr double ABORT_HOVER_ALT = 2; // DECISION-D5 (TRANSPILE-FLAG F01): abort waypoint 2 altitude (m); the lead's spec says 1 m above the pad
    static_assert(ABORT_HOVER_ALT == POLYQP_ABORT_HOVER_ALT_M, "STATIC: the grid capacity derivation (PolyQP_config.h) assumes this altitude");
    // FIX-F02: pin knot-1 accel/jerk to hover thrust along the current b3 (DECISION-D4, off by default)

    // Persist abort calculation state across function calls
    // FIX-F08: 0 = OK, 1 = INFEASIBLE_FLOWN (D2), 2 = FALLBACK (D1)
    // FIX-F08: 0 = OK, 1 = INFEASIBLE_FLOWN (D2), 2 = FALLBACK (D1)
    if (ABORT > 0 && t > ABORT && !state.abortCalculated) {
        // ASTRA-ONLY: the MATLAB Opts (R = 5, InitOxMass/InitFuMass from X_est(14:15), FIX-F29) is not carried: r is hardcoded and ASTRA carries no propellant.
        const TOADWaypointT dummyWP = TOADWaypoint(Vector3d(0, 0, 0));
        TOADWaypointArray &Waypoints = workspace().AbortGuidance.Waypoints; // STATIC: shared workspace
        bool planOk = repmat(dummyWP, 1, 3, Waypoints);

        // FIX-F03: actual thrust axis b3 = C_IB * e3 (quatRot returns C_BI) and db3/dt = w_I x b3 with w_I = C_IB * w_B
        const Vector4d q_est = X_est(R1(1, 4)) / mnorm(Vector4d(X_est(R1(1, 4))));
        const Matrix3_3d C_IB = quatRot(q_est).transpose();
        const Vector3d b3_est = C_IB * Vector3d(0, 0, 1);
        const Vector3d b3dot_est = (C_IB * Vector3d(X_est(R1(11, 13)))).cross(b3_est);
        TOADWaypointOptions options1;
        options1.VelBC = X_est(R1(8, 10));
        options1.VelFixed.setConstant(true);
        options1.Label = "ABT START";
        options1.B3BC = b3_est;
        options1.B3Fixed.setConstant(true);
        options1.B3DotBC = b3dot_est;
        options1.B3DotFixed.setConstant(true);
        Waypoints(1) = TOADWaypoint(X_est(R1(5, 7)), options1); // FIX-F03: knot-1 attitude
        // TRANSPILE-FLAG F02: by default only position and velocity are pinned at knot 1; accel/jerk/snap are free, so thrust and attitude feedforward can step at the abort instant (Phase 1.4). ABORT_PIN_THRUST (DECISION-D4) pins accel and jerk.
        if (ABORT_PIN_THRUST) { // FIX-F02: a0 = (T_hover / m) * b3 - g * e3 and j0 = (T_hover / m) * db3/dt, T_hover / m = g
            options1.AccelBC = ASTRA::g * b3_est - Vector3d(0, 0, ASTRA::g);
            options1.AccelFixed.setConstant(true);
            options1.JerkBC = ASTRA::g * b3dot_est;
            options1.JerkFixed.setConstant(true);
            Waypoints(1) = TOADWaypoint(X_est(R1(5, 7)), options1);
        }
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
        Waypoints(3) = TOADWaypoint(Vector3d(0, 0, 0), options3, Vector3b(false, false, true)); // NO-NAN: [nan; nan; 0], the lateral landing position is free
        // NO-NAN: the MATLAB "non-finite start state -> fallback" branch (FIX-F08) was removed: the estimate is always real
        // TRANSPILE-FLAG F08 (fixed by FIX-F08): MATLAB ignored Report.Feasible and PolyQP errors and had no fallback.
        // TRANSPILE-FLAG F16: the PolyQPStatus (state.planStatus) and Report are returned to the C++ caller; MATLAB keeps only abortStatus. Transpiled as-is per MATLAB.
        // TRANSPILE-FLAG F32 (fixed by FIX-F32): a PolyQP error used to return before abortCalculated was set, so the abort was re-attempted on every later call.
        state.planStatus = PolyQPStatus::Ok;
        if (planOk) {
            InitialGuessT InitialGuess;
            QPTrajT &QPTraj = workspace().AbortGuidance.QPTraj; // STATIC: shared workspace
            // FIX-F08: code generation has no try/catch, so a PolyQP error is caught only in interpreted MATLAB
            // TRANSPILE-FLAG F35: the C++ catches every PolyQP error (including CapacityExceeded) through the status return, like interpreted MATLAB's try/catch.
            // STREAM: Time = InitialGuess.Time(:), States = InitialGuess.X' and Inputs = [U'; U(:, end)'] go straight to Sink row by row
            // STREAM: (sampleOutputs); a failed plan writes nothing, so the fallback below writes the whole trajectory.
            state.planStatus = PolyQP(Waypoints, Sink, InitialGuess, QPTraj, Report);
            planOk = state.planStatus == PolyQPStatus::Ok;
            if (planOk) {
                // FIX-F08: an infeasible plan is flown anyway (D2) and flagged
                if (Report.Feasible) {
                    state.abortStatus = AbortStatus::OK;
                } else {
                    state.abortStatus = AbortStatus::INFEASIBLE_FLOWN;
                }
            }
        }
        if (!planOk) { // FIX-F08: closed-form vertical descent that cannot fail (D1)
            abortFallback(X_est, V1(Waypoints(3).MaxVel, 3), Sink);
            state.abortStatus = AbortStatus::FALLBACK;
        }

        // Turn off abort to avoid recalculating and reference time from
        // when the abort was called
        state.abortCalculated = true; // FIX-F32: set on every path (plan, infeasible plan, fallback): one attempt only
        state.tAbort = t;
    }
    return state.abortStatus;
}

}
