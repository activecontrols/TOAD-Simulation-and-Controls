// Transpiled from TOAD/Guidance/AbortGuidance.m (the abort block formerly in TrajectoryManager.m) @ 5ac6bd23, with the abort fixes of ABORT_FIXES.md
#include "AbortGuidance.h"
#include "PolyQP.h"

namespace PolyQP {

AbortStatus AbortGuidance(Real t, const AbortXEstT &X_est, AbortGuidanceState &state, AbortTimeT &Time, AbortStatesT &States, AbortInputsT &Inputs,
                          ReportT &Report, Real ABORT, bool ABORT_PIN_THRUST) {
    // TRANSPILE-FLAG F09: MATLAB hard-codes the trigger ABORT = 0 (never fires). Transpiled as-is per MATLAB.
    // TRANSPILE-FLAG F31: ABORT is a parameter defaulting to 0 (the MATLAB value) so flight code cannot abort unless a caller sets it; tests set it to 1.
    // TRANSPILE-FLAG F31 (cont.): the persistent tAbort/abortCalculated/abortStatus live in AbortGuidanceState (tAbort = 0, abortCalculated = false, abortStatus = OK). Transpiled as-is per MATLAB.
    // TRANSPILE-FLAG F37: ABORT_PIN_THRUST is a parameter defaulting to the MATLAB constant (false), like ABORT (F31), so tests can switch it on.
    const Real ABORT_HOVER_ALT = 2; // DECISION-D5 (TRANSPILE-FLAG F01): abort waypoint 2 altitude (m); the lead's spec says 1 m above the pad
    // FIX-F02: pin knot-1 accel/jerk to hover thrust along the current b3 (DECISION-D4, off by default)

    // Persist abort calculation state across function calls
    // FIX-F08: 0 = OK, 1 = INFEASIBLE_FLOWN (D2), 2 = FALLBACK (D1)
    // FIX-F08: 0 = OK, 1 = INFEASIBLE_FLOWN (D2), 2 = FALLBACK (D1)
    if (ABORT > 0 && t > ABORT && !state.abortCalculated) {
        // ASTRA-ONLY: the MATLAB Opts (R = 5, InitOxMass/InitFuMass from X_est(14:15), FIX-F29) is not carried: r is hardcoded and ASTRA carries no propellant.
        const TOADWaypointT dummyWP = TOADWaypoint(Vec3(0, 0, 0));
        static TOADWaypointArray Waypoints;
        bool planOk = repmat(dummyWP, 1, 3, Waypoints);

        // FIX-F03: actual thrust axis b3 = C_IB * e3 (quatRot returns C_BI) and db3/dt = w_I x b3 with w_I = C_IB * w_B
        const Vec4 q_est = X_est(R1(1, 4)) / mnorm(Vec4(X_est(R1(1, 4))));
        const Mat3 C_IB = quatRot(q_est).transpose();
        const Vec3 b3_est = C_IB * Vec3(0, 0, 1);
        const Vec3 b3dot_est = (C_IB * Vec3(X_est(R1(11, 13)))).cross(b3_est);
        TOADWaypointOptions options1;
        options1.VelBC = X_est(R1(8, 10));
        options1.Label = "ABT START";
        options1.B3BC = b3_est;
        options1.B3DotBC = b3dot_est;
        Waypoints(1) = TOADWaypoint(X_est(R1(5, 7)), options1); // FIX-F03: knot-1 attitude
        // TRANSPILE-FLAG F02: by default only position and velocity are pinned at knot 1; accel/jerk/snap are free, so thrust and attitude feedforward can step at the abort instant (Phase 1.4). ABORT_PIN_THRUST (DECISION-D4) pins accel and jerk.
        if (ABORT_PIN_THRUST) { // FIX-F02: a0 = (T_hover / m) * b3 - g * e3 and j0 = (T_hover / m) * db3/dt, T_hover / m = g
            options1.AccelBC = ASTRA::g * b3_est - Vec3(0, 0, ASTRA::g);
            options1.JerkBC = ASTRA::g * b3dot_est;
            Waypoints(1) = TOADWaypoint(X_est(R1(5, 7)), options1);
        }
        // TRANSPILE-FLAG F01: abort waypoint 2 altitude is 2 m (ABORT_HOVER_ALT, DECISION-D5); the lead's spec says 1 m above the pad. Transpiled as-is per MATLAB.
        TOADWaypointOptions options2;
        options2.Label = "ABT GND";
        options2.HoldTime = Real(0.0);
        options2.VelBC = Vec3(0, 0, 0);
        options2.AccelBC = Vec3(0, 0, 0);
        options2.JerkBC = Vec3(0, 0, 0);
        Waypoints(2) = TOADWaypoint(Vec3(NaN(), NaN(), ABORT_HOVER_ALT), options2);
        TOADWaypointOptions options3;
        options3.Label = "ABT END";
        options3.HoldTime = Real(0.0);
        options3.MaxVel = Vec3(3, 3, 1);
        options3.StraightLine = true;
        options3.VelBC = Vec3(0, 0, 0);
        options3.AccelBC = Vec3(0, 0, 0);
        options3.JerkBC = Vec3(0, 0, 0);
        Waypoints(3) = TOADWaypoint(Vec3(NaN(), NaN(), 0), options3);
        // FIX-F08: a non-finite start state cannot be planned; it goes straight to the fallback
        planOk = planOk && X_est(R1(5, 10)).allFinite();
        // TRANSPILE-FLAG F08 (fixed by FIX-F08): MATLAB ignored Report.Feasible and PolyQP errors and had no fallback.
        // TRANSPILE-FLAG F16: the PolyQPStatus (state.planStatus) and Report are returned to the C++ caller; MATLAB keeps only abortStatus. Transpiled as-is per MATLAB.
        // TRANSPILE-FLAG F32 (fixed by FIX-F32): a PolyQP error used to return before abortCalculated was set, so the abort was re-attempted on every later call.
        state.planStatus = PolyQPStatus::Ok;
        if (planOk) {
            static InitialGuessT InitialGuess;
            static QPTrajT QPTraj;
            // FIX-F08: code generation has no try/catch, so a PolyQP error is caught only in interpreted MATLAB
            // TRANSPILE-FLAG F35: the C++ catches every PolyQP error (including CapacityExceeded) through the status return, like interpreted MATLAB's try/catch.
            state.planStatus = PolyQP(Waypoints, InitialGuess, QPTraj, Report);
            planOk = state.planStatus == PolyQPStatus::Ok;
            if (planOk) {
                Time = InitialGuess.Time.transpose();
                States = InitialGuess.X.transpose();
                Inputs.resize(InitialGuess.U.cols() + 1, 4);
                Inputs << InitialGuess.U.transpose(), Col1(InitialGuess.U, static_cast<int>(InitialGuess.U.cols())).transpose();
                // FIX-F08: an infeasible plan is flown anyway (D2) and flagged
                if (Report.Feasible) {
                    state.abortStatus = AbortStatus::OK;
                } else {
                    state.abortStatus = AbortStatus::INFEASIBLE_FLOWN;
                }
            }
        }
        if (!planOk) { // FIX-F08: closed-form vertical descent that cannot fail (D1)
            abortFallback(X_est, V1(Waypoints(3).MaxVel, 3), Time, States, Inputs);
            state.abortStatus = AbortStatus::FALLBACK;
        }

        // Turn off abort to avoid recalculating and reference time from
        // when the abort was called
        state.abortCalculated = true; // FIX-F32: set on every path (plan, infeasible plan, fallback): one attempt only
        state.tAbort = t;
    }
    return state.abortStatus;
}

AbortStatus AbortGuidance(float t, const Vector15 &X_est, AbortGuidanceState &state, AbortTimeF &Time, AbortStatesF &States, AbortInputsF &Inputs,
                          ReportT &Report, float ABORT, bool ABORT_PIN_THRUST) {
    const AbortXEstT X_est_real = X_est.cast<Real>();
    static AbortTimeT Time_real;
    static AbortStatesT States_real;
    static AbortInputsT Inputs_real;
    const bool abortCalculatedBefore = state.abortCalculated;
    const AbortStatus st = AbortGuidance(static_cast<Real>(t), X_est_real, state, Time_real, States_real, Inputs_real, Report, static_cast<Real>(ABORT), ABORT_PIN_THRUST);
    if (state.abortCalculated && !abortCalculatedBefore) {
        Time = Time_real.cast<float>();
        States = States_real.cast<float>();
        Inputs = Inputs_real.cast<float>();
    }
    return st;
}

}
