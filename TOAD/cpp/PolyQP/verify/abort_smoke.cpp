// Verification only (never part of the flight build): host smoke test for AbortGuidance.
// ASTRA only: the constants are compiled in (ASTRA_PolyQP_Constants.h).
#include "AbortGuidance.h"
#include "PolyQP.h"

#include <cstdio>

using namespace PolyQP;

static Vector15d makeState(double x, double y, double z, double vx, double vy, double vz) {
    Vector15d X;
    X << 1, 0, 0, 0, x, y, z, vx, vy, vz, 0, 0, 0, 0, 0;
    return X;
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

// STREAM: the smoke test's sink: plain static arrays shaped like the flight trajectory arrays (double here).
static double Time[POLYQP_MAX_NODES];
static double States[POLYQP_MAX_NODES][15];
static double Inputs[POLYQP_MAX_NODES][4];
static int nRows = 0;
static void putRow(void *, int k, double t, const Vector15d &X, const Vector4d &U) {
    Time[k - 1] = t;
    for (int i = 0; i < 15; ++i) States[k - 1][i] = X(i);
    for (int i = 0; i < 4; ++i) Inputs[k - 1][i] = U(i);
    if (k > nRows) nRows = k;
}
static const NodeSink sink{nullptr, POLYQP_MAX_NODES, &putRow};

static void printRow(const char *label, int row) {
    std::printf("  %s:", label);
    for (int j = 0; j < 15; ++j) std::printf(" %.6g", States[row][j]);
    std::printf("\n");
}

int main() {
    int failures = 0;

    {
        AbortGuidanceState state;
        ReportT Report{};
        nRows = 0;
        const AbortStatus st = AbortGuidance(100, makeState(3, -2, 10, 1, 0.5, -1), state, sink, Report);
        const bool ok = st == AbortStatus::OK && !state.abortCalculated && state.tAbort == 0 && nRows == 0;
        std::printf("default ABORT = 0: abortStatus=%s abortCalculated=%d outputs untouched=%d -> %s\n", abortName(st), state.abortCalculated ? 1 : 0, nRows == 0 ? 1 : 0,
                    ok ? "PASS" : "FAIL");
        failures += ok ? 0 : 1;
    }

    const double cases[3][6] = {{3, -2, 10, 1, 0.5, -1}, {0, 0, 5, 0, 0, 0}, {-15, 12, 30, -3, 4, 2}};
    for (int c = 0; c < 3; ++c) {
        AbortGuidanceState state;
        ReportT Report{};
        nRows = 0;
        const Vector15d X_est = makeState(cases[c][0], cases[c][1], cases[c][2], cases[c][3], cases[c][4], cases[c][5]);
        const double ABORT = 1;
        const double t = 2;
        const AbortStatus ast = AbortGuidance(t, X_est, state, sink, Report, ABORT);
        const PolyQPStatus st = state.planStatus;
        std::printf("case %d: X_est pos=(%g, %g, %g) vel=(%g, %g, %g)\n", c + 1, cases[c][0], cases[c][1], cases[c][2], cases[c][3], cases[c][4], cases[c][5]);
        std::printf("  abortStatus=%s planStatus=%s abortCalculated=%d tAbort=%g Feasible=%d T_total=%.9g samples=%d OuterIterations=%d\n", abortName(ast), statusName(st),
                    state.abortCalculated ? 1 : 0, state.tAbort, Report.Feasible ? 1 : 0, Report.T_total, nRows, Report.OuterIterations);
        if (st != PolyQPStatus::Ok || ast == AbortStatus::FALLBACK) {
#if defined(POLYQP_ERROR_MESSAGES)
            std::printf("  error %s: %s\n", lastError().identifier, lastError().message);
#endif
            ++failures;
            continue;
        }
        printRow("States(1,:)  ", 0);
        printRow("States(end,:)", nRows - 1);
        bool dup = true;
        for (int i = 0; i < 4; ++i) dup = dup && Inputs[nRows - 1][i] == Inputs[nRows - 2][i];
        const bool timeOk = Time[0] == 0 && Time[nRows - 1] == Report.T_total;
        const int n0 = nRows;
        const AbortStatus st2 = AbortGuidance(3, X_est, state, sink, Report, ABORT);
        const bool once = st2 == ast && state.tAbort == t && nRows == n0;
        std::printf("  last-input-duplicated=%d time-grid=%d no-recompute=%d\n", dup ? 1 : 0, timeOk ? 1 : 0, once ? 1 : 0);
        failures += (dup && timeOk && once) ? 0 : 1;
    }

    // FIX-F08 fallback through CapacityExceeded: a 300 m start needs more output nodes than POLYQP_MAX_NODES.
    {
        AbortGuidanceState state;
        ReportT Report{};
        nRows = 0;
        const Vector15d X_est = makeState(3, -2, 300, 1, 0.5, -1);
        const AbortStatus st = AbortGuidance(2, X_est, state, sink, Report, 1);
        const int n = nRows;
        const bool start = n >= 2 && States[0][4] == X_est(4) && States[0][5] == X_est(5) && States[0][6] == X_est(6) && Time[0] == 0;
        const bool end = n >= 2 && States[n - 1][6] == 0 && States[n - 1][9] == 0;
        const bool ok = st == AbortStatus::FALLBACK && state.planStatus == PolyQPStatus::CapacityExceeded && state.abortCalculated && start && end && n <= POLYQP_MAX_NODES;
        std::printf("fallback (300 m start, CapacityExceeded): abortStatus=%s planStatus=%s samples=%d T=%.6g start=%d end-at-rest=%d -> %s\n", abortName(st),
                    statusName(state.planStatus), n, n >= 1 ? Time[n - 1] : -1.0, start ? 1 : 0, end ? 1 : 0, ok ? "PASS" : "FAIL");
        failures += ok ? 0 : 1;
    }

    std::printf("%s (%d failure(s))\n", failures == 0 ? "SMOKE PASS" : "SMOKE FAIL", failures);
    return failures == 0 ? 0 : 1;
}
