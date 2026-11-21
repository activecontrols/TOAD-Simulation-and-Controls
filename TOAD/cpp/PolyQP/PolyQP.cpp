// Transpiled from TOAD/Guidance/PolynomialQP/PolyQP.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

static PolyQPStatus solveAllLegs(PlanT &Plan, const Matrix10_10d &A0, const Matrix10_10d &H0, const LimitsT &Limits, bool doAllocateTimes, QPTrajT &QPTraj);
static PolyQPStatus buildReport(bool Feasible, const Vector9d &MaxUtilPerType, int outerIter, const QPTrajT &QPTraj, ReportT &Report);

// STATIC: the one shared workspace of the solver (see PolyQP_types.h); not reentrant.
static PolyQPWorkspace g_workspace;
PolyQPWorkspace &workspace() { return g_workspace; }

PolyQPStatus PolyQP(const TOADWaypointArray &Waypoints, const NodeSink &Sink, InitialGuessT &InitialGuess, QPTrajT &QPTraj, ReportT &Report) {
    // Waypoint-Driven Closed-Form Polynomial Trajectory Generator
    // TRANSPILE-FLAG F10: the nargin < 2 branch (Constants = []) cannot occur in C++ (the ASTRA constants are compiled in) and is dead in MATLAB (buildVehicleLimits dereferences Constants.Vehicle). Transpiled as-is per MATLAB.
    // TRANSPILE-FLAG F33: the MATLAB Opts argument is not carried: r is the constexpr POLYQP_R = 5 (the MATLAB default), and InitOxMass/InitFuMass were superseded (ASTRA carries no propellant).
    const int r = POLYQP_R;

    // 1. Build vehicle limits and resolve waypoints
    LimitsT Limits;
    POLYQP_TRY(buildVehicleLimits(Limits));
    PlanT &Plan = workspace().PolyQP.Plan; // STATIC: shared workspace
    POLYQP_TRY(resolveWaypoints(Waypoints, Limits, Plan));
    Matrix10_10d A0;
    Matrix10_10d Q0Ignored;
    Matrix10_10d H0;
    POLYQP_TRY(assembleQP(A0, Q0Ignored, H0));

    // 2. Convert CG position BCs to Centre-of-Percussion P-space and solve all legs
    // Run Richter time redistribution only on legs where intermediate knots have no fixed VelBC
    POLYQP_TRY(solveAllLegs(Plan, A0, H0, Limits, true, QPTraj));

    // 3. If no vehicle constants provided, return trajectory only
    // TRANSPILE-FLAG F10: Limits.HasVehicle is always true (buildVehicleLimits dereferences Constants.Vehicle), so this branch is dead. Transpiled as-is per MATLAB.
    if (!Limits.HasVehicle) {
        InitialGuess.Nn = 0;
        InitialGuess.T_total = 0;
        InitialGuess.dt = 0;
        POLYQP_TRY(buildReport(true, Vector9d::Zero(), 0, QPTraj, Report));
        return PolyQPStatus::Ok;
    }

    // 4. Unified Evaluate -> Bidirectional Leg Time Scaling Loop
    const int MaxOuter = 40;
    int stallCount = 0;
    int outerIter = 0;

    const int nLegs = Plan.LegCount;
    const int nSeg = Plan.SegCount; // STATIC: number of valid columns of EvalOut.SegUtil
    Array1_Sb legCanContract = Array1_Sb::Constant(true);
    for (int leg = 1; leg <= nLegs; ++leg) {
        const int m1 = V1(Plan.LegStartSeg, leg);
        const int m2 = V1(Plan.LegEndSeg, leg);
        // Only contract legs whose fixed derivative BCs (order >= 1) are homogeneous (0 or NaN)
        Array3_5_Kd &dLeg = workspace().PolyQP.dLeg;
        sub3(Plan.FixCg, 1, 3, 2, r, m1, m2 + 1, dLeg);
        Array3_5_Kb &dLegFixed = workspace().PolyQP.dLegFixed;
        sub3(Plan.FixCgFixed, 1, 3, 2, r, m1, m2 + 1, dLegFixed);
        // TRANSPILE-FLAG F05: false for every abort (knot 1 carries the nonzero estimated velocity), so the abort never contracts below its initial time guess. Transpiled as-is per MATLAB.
        bool anyNonzeroFixed = false; // NO-NAN: any(dLeg(:) ~= 0 & ~isnan(dLeg(:)))
        for (int q = 0; q < dLeg.numel(); ++q) {
            if (dLegFixed.data[q] && dLeg.data[q] != double(0)) anyNonzeroFixed = true;
        }
        if (anyNonzeroFixed) {
            V1(legCanContract, leg) = false;
        }
    }
    Array1_Sb legHasDilated = Array1_Sb::Constant(false);

    EvalOutT &EvalOut = workspace().PolyQP.EvalOut;
    bool Feasible;
    POLYQP_TRY(evaluateConstraints(QPTraj, Plan, Limits, EvalOut, Feasible));
    double bestUtil = EvalOut.MaxOverallUtil;
    Matrix1_Sd bestSegT = Plan.SegT;
    Matrix1_Sd prevSegUtil = Matrix1_Sd::Zero();
    maxDim1(EvalOut.SegUtil.leftCols(nSeg), prevSegUtil);

    while (outerIter < MaxOuter) {
        bool didScale = false;
        bool didContract = false;

        for (int leg = 1; leg <= nLegs; ++leg) {
            const int m1 = V1(Plan.LegStartSeg, leg);
            const int m2 = V1(Plan.LegEndSeg, leg);
            const double s_leg = maxv(EvalOut.SegDilate(R1(m1, m2)));
            if (s_leg > double(1.0)) {
                const double s_apply = std::min(double(2.0), double(1.02) * s_leg);
                Plan.SegT(R1(m1, m2)) = Plan.SegT(R1(m1, m2)) * s_apply;
                V1(legHasDilated, leg) = true;
                didScale = true;
            } else if (Feasible && V1(legCanContract, leg) && !V1(legHasDilated, leg) && (double(1.02) * s_leg < double(0.97))) {
                const double s_min_leg = maxv((Limits.MinSegTime / Plan.SegT(R1(m1, m2)).array()).matrix());
                const double s_apply = maxv({double(0.60), s_min_leg, double(1.02) * s_leg});
                if (s_apply < double(1.0) - double(1e-4)) {
                    Plan.SegT(R1(m1, m2)) = Plan.SegT(R1(m1, m2)) * s_apply;
                    didScale = true;
                    didContract = true;
                }
            }
        }

        if (!didScale) {
            break;
        }

        outerIter = outerIter + 1;

        // Re-solve legs with updated segment durations (preserving optimal relative time shares)
        Plan.R = r;
        POLYQP_TRY(solveAllLegs(Plan, A0, H0, Limits, false, QPTraj));
        POLYQP_TRY(evaluateConstraints(QPTraj, Plan, Limits, EvalOut, Feasible));

        if (didContract) {
            bestUtil = EvalOut.MaxOverallUtil;
            if (Feasible) {
                bestSegT = Plan.SegT;
            }
            maxDim1(EvalOut.SegUtil.leftCols(nSeg), prevSegUtil);
            stallCount = 0;
            continue;
        }

        // Stall detection during dilation phase
        Matrix1_Sd curSegUtil = Matrix1_Sd::Zero();
        maxDim1(EvalOut.SegUtil.leftCols(nSeg), curSegUtil);
        const Array1_Sb violMask = curSegUtil.array() > double(1.0) + double(1e-3);
        bool anyViolImproved = false;
        for (int q = 1; q <= nSeg; ++q) {
            if (V1(violMask, q) && V1(curSegUtil, q) < V1(prevSegUtil, q) - double(1e-4)) anyViolImproved = true;
        }
        const bool madeProgress = (EvalOut.MaxOverallUtil < bestUtil - double(1e-4)) || anyViolImproved;
        prevSegUtil = curSegUtil;

        if (EvalOut.MaxOverallUtil < bestUtil - double(1e-4)) {
            bestUtil = EvalOut.MaxOverallUtil;
            bestSegT = Plan.SegT;
        }

        if (madeProgress) {
            stallCount = 0;
        } else {
            stallCount = stallCount + 1;
            if (stallCount >= 4) {
                Plan.SegT = bestSegT;
                POLYQP_TRY(solveAllLegs(Plan, A0, H0, Limits, false, QPTraj));
                POLYQP_TRY(evaluateConstraints(QPTraj, Plan, Limits, EvalOut, Feasible));
                break;
            }
        }
    }

    // 5. Sample uniform output grid and build Report
    POLYQP_TRY(sampleOutputs(QPTraj, Plan, Limits, Sink, InitialGuess));
    POLYQP_TRY(buildReport(Feasible, EvalOut.MaxUtilPerType, outerIter, QPTraj, Report));
    return PolyQPStatus::Ok;
}

static PolyQPStatus solveAllLegs(PlanT &Plan, const Matrix10_10d &A0, const Matrix10_10d &H0, const LimitsT &Limits, bool doAllocateTimes, QPTrajT &QPTraj) {
    // Solve All Flight Legs with P-Space Fixed-Point Refinement
    const int MAX_SEG = POLYQP_MAX_SEG;
    const int MAX_KNOTS = MAX_SEG + 1;
    const int r = Plan.R;
    const int nc = 2 * r;
    const int M = Plan.SegCount;
    const int K = Plan.KnotCount;
    const double rho = Limits.Rho;
    const double g = Limits.g;

    // Initialize P-space boundary conditions: P = r_cg + rho * b3 (init b3 = e3)
    // TRANSPILE-FLAG F03 (fixed by FIX-F03): the knot-1 P-offset assumed b3 = e3 (rho ~ 0.22 m for ASTRA) and the velocity BC was applied to dP/dt, not v_cg + rho*db3/dt; a knot with B3BC/B3DotBC now uses the given attitude.
    PolyQPWorkspace::SolveAllLegsT &W = workspace().solveAllLegs; // STATIC: shared workspace
    Array3_5_Kd &FixP = W.FixP;
    sub3(Plan.FixCg, 1, 3, 1, r, 1, K, FixP);
    Array3_5_Kb &FixPFixed = W.FixPFixed; // NO-NAN: ~isnan(FixP)
    sub3(Plan.FixCgFixed, 1, 3, 1, r, 1, K, FixPFixed);
    for (int k = 1; k <= K; ++k) {
        const bool hasB3 = Col1(Plan.B3Fixed, k).all(); // FIX-F03: a knot with a given attitude uses P = r_cg + rho * b3
        const bool hasB3Dot = Col1(Plan.B3DotFixed, k).all(); // FIX-F03: and dP/dt = v_cg + rho * db3/dt where the velocity is fixed
        for (int d = 1; d <= 3; ++d) {
            if (FixPFixed(d, 1, k)) {
                if (hasB3) { // FIX-F03: P = r_cg + rho * b3 with the given attitude
                    FixP(d, 1, k) = Plan.FixCg(d, 1, k) + rho * M1(Plan.B3BC, d, k);
                } else if (d == 3) {
                    FixP(d, 1, k) = Plan.FixCg(d, 1, k) + rho;
                } else {
                    FixP(d, 1, k) = Plan.FixCg(d, 1, k);
                }
            }
            if (hasB3Dot && FixPFixed(d, 2, k)) { // FIX-F03: dP/dt = v_cg + rho * db3/dt
                FixP(d, 2, k) = Plan.FixCg(d, 2, k) + rho * M1(Plan.B3DotBC, d, k);
            }
        }
    }

    Array10_S_3d &CoefAll = W.CoefAll;
    CoefAll.setConstant(POLYQP_MAX_NC, MAX_SEG, 3, double(0));
    Array3_5_Kd &DerivAll = W.DerivAll;
    DerivAll = FixP;
    Array3_5_Kb &DerivAllFixed = W.DerivAllFixed;
    DerivAllFixed = FixPFixed;
    double totalCost = 0;

    // Optional one-time Richter time redistribution on multi-segment legs where intermediate knots have at least one free velocity axis
    if (doAllocateTimes) {
        for (int leg = 1; leg <= Plan.LegCount; ++leg) {
            const int m1 = V1(Plan.LegStartSeg, leg);
            const int m2 = V1(Plan.LegEndSeg, leg);
            if (m2 > m1) {
                // NO-NAN: all(any(isnan(midVelBins), 1)) with midVelBins = FixP(:, 2, (m1 + 1):m2): every intermediate knot has a free velocity axis
                bool allAnyFree = true;
                for (int k = m1 + 1; k <= m2; ++k) {
                    if (FixPFixed(1, 2, k) && FixPFixed(2, 2, k) && FixPFixed(3, 2, k)) allAnyFree = false;
                }
                if (allAnyFree) {
                    Array3_5_Kd &FixLeg = W.FixLegAT;
                    sub3(FixP, 1, 3, 1, r, m1, m2 + 1, FixLeg);
                    Array3_5_Kb &FixLegFixed = W.FixLegFixedAT;
                    sub3(FixPFixed, 1, 3, 1, r, m1, m2 + 1, FixLegFixed);
                    Matrix1_Sd T_leg = Matrix1_Sd::Zero();
                    T_leg.head(m2 - m1 + 1) = Plan.SegT(R1(m1, m2));
                    Matrix1_Sd T_opt = Matrix1_Sd::Zero();
                    double CostIgnored;
                    POLYQP_TRY(allocateTimes(FixLeg, FixLegFixed, T_leg, A0, H0, Limits.MinSegTime, T_opt, CostIgnored));
                    Plan.SegT(R1(m1, m2)) = T_opt.head(m2 - m1 + 1);
                }
            }
        }
    }

    // Fixed-point P-space BC refinement at pass-through knots (<= 3 passes)
    const bool hasPassThrough = (!Plan.KnotIsRest(R1(1, K))).any();
    int nPasses = 1;
    if (hasPassThrough && Limits.HasVehicle) {
        // TRANSPILE-FLAG F04: for the abort, passes 2-3 recompute pass 1 exactly (the only pass-through knot is knot 1, which is never refined). Transpiled as-is per MATLAB.
        nPasses = 3;
    }

    for (int pass = 1; pass <= nPasses; ++pass) {
        totalCost = 0;

        // Solve each decoupled flight leg
        for (int leg = 1; leg <= Plan.LegCount; ++leg) {
            const int m1 = V1(Plan.LegStartSeg, leg);
            const int m2 = V1(Plan.LegEndSeg, leg);
            Array3_5_Kd &FixLeg = W.FixLeg;
            sub3(FixP, 1, 3, 1, r, m1, m2 + 1, FixLeg);
            Array3_5_Kb &FixLegFixed = W.FixLegFixed;
            sub3(FixPFixed, 1, 3, 1, r, m1, m2 + 1, FixLegFixed);
            Matrix1_Sd T_leg = Matrix1_Sd::Zero();
            T_leg.head(m2 - m1 + 1) = Plan.SegT(R1(m1, m2));

            Array10_S_3d &CoefLeg = W.CoefLeg;
            Array3_5_Kd &DerivLeg = W.DerivLeg;
            double costLeg;
            POLYQP_TRY(solveLeg(FixLeg, FixLegFixed, T_leg, A0, H0, CoefLeg, DerivLeg, costLeg));
            setSub3(CoefAll, 1, nc, m1, m2, 1, 3, CoefLeg);
            setSub3(DerivAll, 1, 3, 1, r, m1, m2 + 1, DerivLeg);
            fillSub3(DerivAllFixed, 1, 3, 1, r, m1, m2 + 1, true); // NO-NAN: solveLeg returns every derivative of the leg
            totalCost = totalCost + costLeg;
        }

        // Populate Hold and Cruise segment polynomial coefficients
        for (int m = 1; m <= M; ++m) {
            if (V1(Plan.IsHold, m)) {
                for (int d = 1; d <= 3; ++d) {
                    fillSub3(CoefAll, 1, nc, m, m, d, d, double(0));
                    CoefAll(1, m, d) = FixP(d, 1, m);
                }
            } else if (V1(Plan.IsCruise, m)) {
                for (int d = 1; d <= 3; ++d) {
                    fillSub3(CoefAll, 1, nc, m, m, d, d, double(0));
                    CoefAll(1, m, d) = FixP(d, 1, m);
                    CoefAll(2, m, d) = FixP(d, 2, m);
                }
            }
        }

        // Update P = r_cg + rho * b3 at pass-through knots using solved acceleration
        if (pass < nPasses) {
            // TRANSPILE-FLAG F03 (fixed by FIX-F03): refinement covers k = 2:(K-1) only; the knot-1 attitude now comes from B3BC instead of the b3 = e3 assumption.
            for (int k = 2; k <= K - 1; ++k) {
                if (!V1(Plan.KnotIsRest, k) && !Col1(Plan.B3Fixed, k).all()) { // FIX-F03: a given knot attitude is not overwritten by the refinement
                    const Vector3d a_k = getFiber3(DerivAll, 3, k);
                    const Vector3d e_k = a_k + Vector3d(0, 0, g);
                    const double ne = mnorm(e_k);
                    if (ne > double(1e-6)) {
                        const Vector3d b3_k = e_k / ne;
                        for (int d = 1; d <= 3; ++d) {
                            if (Plan.FixCgFixed(d, 1, k)) {
                                FixP(d, 1, k) = Plan.FixCg(d, 1, k) + rho * V1(b3_k, d);
                            }
                        }
                    }
                }
            }
        }
    }

    Matrix1_Kd KeyTimes = Matrix1_Kd::Zero();
    Matrix1_Sd cumsumSegT = Matrix1_Sd::Zero();
    cumsum(Plan.SegT(R1(1, M)), cumsumSegT);
    V1(KeyTimes, 1) = double(0);
    KeyTimes(R1(2, M + 1)) = cumsumSegT.head(M);

    Array3_6_Kd &FixOut = W.FixOut;
    FixOut.setConstant(3, 6, MAX_KNOTS, double(0));
    setSub3(FixOut, 1, 3, 1, r, 1, K, DerivAll);
    Array3_6_Kb &FixOutFixed = W.FixOutFixed; // NO-NAN: FixOut = nan(3, 6, MAX_KNOTS) before the solved derivatives are written
    FixOutFixed.setConstant(3, 6, MAX_KNOTS, false);
    setSub3(FixOutFixed, 1, 3, 1, r, 1, K, DerivAllFixed);

    QPTraj.R = r;
    QPTraj.Count = M;
    QPTraj.T = Plan.SegT;
    QPTraj.KeyTimes = KeyTimes;
    sub3(CoefAll, 1, nc, 1, CoefAll.size(2), 1, CoefAll.size(3), QPTraj.Coef);
    QPTraj.IsHold = Plan.IsHold;
    QPTraj.IsCruise = Plan.IsCruise;
    QPTraj.Fix = FixOut;
    QPTraj.FixFixed = FixOutFixed;
    QPTraj.Cost = totalCost;
    QPTraj.T_total = sumv(Plan.SegT(R1(1, M)));
    QPTraj.FlatOffset = rho;
    return PolyQPStatus::Ok;
}

static PolyQPStatus buildReport(bool Feasible, const Vector9d &MaxUtilPerType, int outerIter, const QPTrajT &QPTraj, ReportT &Report) {
    // Build Diagnostic Output Report
    const int M = QPTraj.Count;
    Report.Feasible = Feasible;
    Report.MaxOverallUtil = maxv(MaxUtilPerType);
    Report.UtilVel = V1(MaxUtilPerType, 1);
    Report.UtilAccel = V1(MaxUtilPerType, 2);
    Report.UtilThrustHi = V1(MaxUtilPerType, 3);
    Report.UtilThrustLo = V1(MaxUtilPerType, 4);
    Report.UtilGimbal = V1(MaxUtilPerType, 5);
    Report.UtilGimbalRate = V1(MaxUtilPerType, 6);
    Report.UtilThrustRate = V1(MaxUtilPerType, 7);
    Report.UtilRoll = V1(MaxUtilPerType, 8);
    Report.UtilTilt = V1(MaxUtilPerType, 9);
    Report.OuterIterations = outerIter;
    Report.Cost = QPTraj.Cost;
    Report.SegTimes = Matrix1_Sd::Zero();
    Report.SegTimes.head(M) = QPTraj.T.head(M);
    Report.T_total = QPTraj.T_total;
    return PolyQPStatus::Ok;
}

}
