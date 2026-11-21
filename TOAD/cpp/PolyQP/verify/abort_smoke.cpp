// Verification only (never part of the flight build): host smoke test for AbortGuidance.
// ASTRA constants copied from LoadTOADParams(0) @ 5ac6bd23.
#include "AbortGuidance.h"
#include "PolyQP.h"
#include "matlab_compat.h"

#include <cstdio>
#include <limits>

using namespace PolyQP;

static PolyQPConstants astraConstants() {
    PolyQPConstants c;
    c.Vehicle = 0;
    c.m_dry = 1.1690;
    c.g = 9.80145;
    c.rTB = 0.26;
    c.J = Vec3(0.067, 0.067, 0.02).asDiagonal();
    c.MaxThrust = 16;
    c.MaxMdot = 0;
    c.OF = 1;
    c.OxMass = 0;
    c.FuMass = 0;
    c.OxHeight = 1;
    c.FuHeight = 1;
    c.OxRadius = 1;
    c.FuRadius = 1;
    c.Ox_Z = 1;
    c.Fu_Z = 1;
    return c;
}

static AbortXEstT makeState(Real x, Real y, Real z, Real vx, Real vy, Real vz) {
    AbortXEstT X(15);
    X << 1, 0, 0, 0, x, y, z, vx, vy, vz, 0, 0, 0, 0, 0;
    return X;
}

static void printRow(const char *label, const AbortStatesT &S, Eigen::Index row) {
    std::printf("  %s:", label);
    for (Eigen::Index j = 0; j < S.cols(); ++j) std::printf(" %.6g", static_cast<double>(S(row, j)));
    std::printf("\n");
}

static const char *statusName(PolyQPStatus s) {
    switch (s) {
    case PolyQPStatus::Ok: return "Ok";
    case PolyQPStatus::InvalidWaypointCount: return "InvalidWaypointCount";
    case PolyQPStatus::InvalidOrder: return "InvalidOrder";
    case PolyQPStatus::InvalidHDGRef: return "InvalidHDGRef";
    case PolyQPStatus::InconsistentHDGRef: return "InconsistentHDGRef";
    case PolyQPStatus::InvalidLimits: return "InvalidLimits";
    case PolyQPStatus::IllPosedQP: return "IllPosedQP";
    case PolyQPStatus::SingularRPP: return "SingularRPP";
    case PolyQPStatus::CapacityExceeded: return "CapacityExceeded";
    }
    return "?";
}

static const char *abortName(AbortStatus s) {
    switch (s) {
    case AbortStatus::OK: return "OK";
    case AbortStatus::INFEASIBLE_FLOWN: return "INFEASIBLE_FLOWN";
    case AbortStatus::FALLBACK: return "FALLBACK";
    }
    return "?";
}

static AbortTimeT Time;
static AbortStatesT States;
static AbortInputsT Inputs;

int main() {
    const PolyQPConstants constantsTOAD = astraConstants();
    int failures = 0;

    {
        AbortGuidanceState state;
        ReportT Report{};
        Time.resize(1);
        Time(0) = -1;
        const AbortStatus st = AbortGuidance(Real(100), makeState(3, -2, 10, 1, 0.5, -1), constantsTOAD, state, Time, States, Inputs, Report);
        const bool ok = st == AbortStatus::OK && !state.abortCalculated && state.tAbort == 0 && Time.size() == 1 && Time(0) == -1;
        std::printf("default ABORT = 0: abortStatus=%s abortCalculated=%d outputs untouched=%d -> %s\n", abortName(st), state.abortCalculated ? 1 : 0,
                    (Time.size() == 1 && Time(0) == -1) ? 1 : 0, ok ? "PASS" : "FAIL");
        failures += ok ? 0 : 1;
    }

    const Real cases[3][6] = {{3, -2, 10, 1, 0.5, -1}, {0, 0, 5, 0, 0, 0}, {-15, 12, 30, -3, 4, 2}};
    for (int c = 0; c < 3; ++c) {
        AbortGuidanceState state;
        ReportT Report{};
        const AbortXEstT X_est = makeState(cases[c][0], cases[c][1], cases[c][2], cases[c][3], cases[c][4], cases[c][5]);
        const Real ABORT = 1;
        const Real t = 2;
        const AbortStatus ast = AbortGuidance(t, X_est, constantsTOAD, state, Time, States, Inputs, Report, ABORT);
        const PolyQPStatus st = state.planStatus;
        std::printf("case %d: X_est pos=(%g, %g, %g) vel=(%g, %g, %g)\n", c + 1, static_cast<double>(cases[c][0]), static_cast<double>(cases[c][1]),
                    static_cast<double>(cases[c][2]), static_cast<double>(cases[c][3]), static_cast<double>(cases[c][4]), static_cast<double>(cases[c][5]));
        std::printf("  abortStatus=%s planStatus=%s abortCalculated=%d tAbort=%g Feasible=%d T_total=%.9g samples=%ld OuterIterations=%d Inputs=%ldx%ld\n", abortName(ast), statusName(st),
                    state.abortCalculated ? 1 : 0, static_cast<double>(state.tAbort), Report.Feasible ? 1 : 0, static_cast<double>(Report.T_total),
                    static_cast<long>(Time.size()), Report.OuterIterations, static_cast<long>(Inputs.rows()), static_cast<long>(Inputs.cols()));
        if (st != PolyQPStatus::Ok || ast == AbortStatus::FALLBACK) {
#if defined(POLYQP_ERROR_MESSAGES)
            std::printf("  error %s: %s\n", lastError().identifier, lastError().message);
#endif
            ++failures;
            continue;
        }
        printRow("States(1,:)  ", States, 0);
        printRow("States(end,:)", States, States.rows() - 1);
        const bool shapes = Time.size() == States.rows() && Inputs.rows() == States.rows() && States.cols() == 15 && Inputs.cols() == 4 &&
                            Inputs.row(Inputs.rows() - 1) == Inputs.row(Inputs.rows() - 2);
        const Real T_last = Time(Time.size() - 1);
        const bool timeOk = Time(0) == 0 && T_last == Report.T_total;
        const Real T0 = Report.T_total;
        const AbortStatus st2 = AbortGuidance(Real(3), X_est, constantsTOAD, state, Time, States, Inputs, Report, ABORT);
        const bool once = st2 == ast && state.tAbort == t && Report.T_total == T0;
        std::printf("  shapes=%d time-grid=%d no-recompute=%d\n", shapes ? 1 : 0, timeOk ? 1 : 0, once ? 1 : 0);
        failures += (shapes && timeOk && once) ? 0 : 1;
    }

    {
        AbortGuidanceState state;
        ReportT Report{};
        Vector15 X_est_f;
        X_est_f << 1.0f, 0.0f, 0.0f, 0.0f, 3.0f, -2.0f, 10.0f, 1.0f, 0.5f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
        static AbortTimeF TimeF;
        static AbortStatesF StatesF;
        static AbortInputsF InputsF;
        const AbortStatus st = AbortGuidance(2.0f, X_est_f, constantsTOAD, state, TimeF, StatesF, InputsF, Report, 1.0f);
        const bool ok = st != AbortStatus::FALLBACK && state.planStatus == PolyQPStatus::Ok && state.abortCalculated && TimeF.size() == StatesF.rows() &&
                        StatesF.cols() == 15;
        std::printf("float flight overload: abortStatus=%s samples=%ld T_total=%.9g -> %s\n", abortName(st), static_cast<long>(TimeF.size()),
                    static_cast<double>(Report.T_total), ok ? "PASS" : "FAIL");
        failures += ok ? 0 : 1;
    }

    // FIX-F08 fallback: a non-finite estimate and a forced PolyQP error (z0 = 1e300 makes T_m^(1-2r) underflow -> SingularRPP)
    for (int c = 0; c < 2; ++c) {
        AbortGuidanceState state;
        ReportT Report{};
        AbortXEstT X_est = makeState(3, -2, 10, 1, 0.5, -1);
        if (c == 0) {
            X_est(8) = std::numeric_limits<Real>::quiet_NaN();
        } else {
            X_est(6) = Real(1e300);
        }
        const AbortStatus st = AbortGuidance(Real(2), X_est, constantsTOAD, state, Time, States, Inputs, Report, Real(1));
        const Eigen::Index n = States.rows();
        const bool shapes = Time.size() == n && Inputs.rows() == n && n >= 2;
        const bool start = shapes && States(0, 4) == X_est(4) && States(0, 5) == X_est(5) && States(0, 6) == X_est(6) && Time(0) == 0;
        const bool end = shapes && States(n - 1, 6) == 0 && States(n - 1, 9) == 0 && States.allFinite() && Inputs.allFinite();
        const Eigen::Index n0 = n;
        const AbortStatus st2 = AbortGuidance(Real(3), X_est, constantsTOAD, state, Time, States, Inputs, Report, Real(1));
        const bool once = st2 == AbortStatus::FALLBACK && States.rows() == n0 && state.tAbort == Real(2);
        const bool ok = st == AbortStatus::FALLBACK && state.abortCalculated && start && end && once;
        std::printf("fallback (%s): abortStatus=%s planStatus=%s samples=%ld T=%.6g start=%d end-at-rest=%d no-recompute=%d -> %s\n",
                    c == 0 ? "NaN velocity" : "forced PolyQP error", abortName(st), statusName(state.planStatus), static_cast<long>(n),
                    shapes ? static_cast<double>(Time(n - 1)) : -1.0, start ? 1 : 0, end ? 1 : 0, once ? 1 : 0, ok ? "PASS" : "FAIL");
        failures += ok ? 0 : 1;
    }

    std::printf("%s (%d failure(s))\n", failures == 0 ? "SMOKE PASS" : "SMOKE FAIL", failures);
    return failures == 0 ? 0 : 1;
}
