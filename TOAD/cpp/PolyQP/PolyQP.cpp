// Transpiled from TOAD/Guidance/PolynomialQP/PolyQP.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

static PolyQPStatus solveAllLegs(PlanT &Plan, const MatNC &A0, const MatNC &H0, const LimitsT &Limits, bool doAllocateTimes, QPTrajT &QPTraj);
static PolyQPStatus buildReport(bool Feasible, const Vec9 &MaxUtilPerType, int outerIter, const QPTrajT &QPTraj, ReportT &Report);

PolyQPStatus PolyQP(const TOADWaypointArray &Waypoints, const PolyQPConstants &Constants, InitialGuessT &InitialGuess, QPTrajT &QPTraj, ReportT &Report,
                    const PolyQPOpts &Opts) {
    // Waypoint-Driven Closed-Form Polynomial Trajectory Generator
    // TRANSPILE-FLAG F10: the nargin < 2 branch (Constants = []) cannot occur in C++ (Constants is always supplied) and is dead in MATLAB (buildVehicleLimits dereferences Constants.Vehicle). Transpiled as-is per MATLAB.
    // TRANSPILE-FLAG F33: nargin < 3 / isempty(Opts) is the default argument Opts = PolyQPOpts() (R = 5, placed last); Opts.R is always present, so ~isfield/isempty(Opts.R) is constant false. Transpiled as-is per MATLAB.
    const bool OptsRMissing = false;
    int r;
    if (OptsRMissing) {
        r = 5;
    } else {
        r = Opts.R;
    }

    // 1. Build vehicle limits and resolve waypoints
    LimitsT Limits;
    POLYQP_TRY(buildVehicleLimits(Constants, Limits, Opts)); // FIX-F29: Opts.InitOxMass/InitFuMass set the initial propellant
    static PlanT Plan;
    POLYQP_TRY(resolveWaypoints(Waypoints, r, Limits, Plan));
    MatNC A0;
    MatNC Q0Ignored;
    MatNC H0;
    POLYQP_TRY(assembleQP(r, A0, Q0Ignored, H0));

    // 2. Convert CG position BCs to Centre-of-Percussion P-space and solve all legs
    // Run Richter time redistribution only on legs where intermediate knots have no fixed VelBC
    POLYQP_TRY(solveAllLegs(Plan, A0, H0, Limits, true, QPTraj));

    // 3. If no vehicle constants provided, return trajectory only
    // TRANSPILE-FLAG F10: Limits.HasVehicle is always true (buildVehicleLimits dereferences Constants.Vehicle), so this branch is dead. Transpiled as-is per MATLAB.
    if (!Limits.HasVehicle) {
        InitialGuess.Time.resize(0);
        InitialGuess.X.resize(15, 0);
        InitialGuess.U.resize(4, 0);
        InitialGuess.T_total = 0;
        InitialGuess.dt_row.resize(0);
        POLYQP_TRY(buildReport(true, Vec9::Zero(9), 0, QPTraj, Report));
        return PolyQPStatus::Ok;
    }

    // 4. Unified Evaluate -> Bidirectional Leg Time Scaling Loop
    const int MaxOuter = 40;
    int stallCount = 0;
    int outerIter = 0;

    const int nLegs = Plan.LegCount;
    BoolSeg legCanContract = BoolSeg::Constant(nLegs, true);
    for (int leg = 1; leg <= nLegs; ++leg) {
        const int m1 = V1(Plan.LegStartSeg, leg);
        const int m2 = V1(Plan.LegEndSeg, leg);
        // Only contract legs whose fixed derivative BCs (order >= 1) are homogeneous (0 or NaN)
        static Fix3 dLeg;
        sub3(Plan.FixCg, 1, 3, 2, r, m1, m2 + 1, dLeg);
        // TRANSPILE-FLAG F05: false for every abort (knot 1 carries the nonzero estimated velocity), so the abort never contracts below its initial time guess. Transpiled as-is per MATLAB.
        if (((colon(dLeg) != Real(0)) && !colon(dLeg).isNaN()).any()) {
            V1(legCanContract, leg) = false;
        }
    }
    BoolSeg legHasDilated = BoolSeg::Constant(nLegs, false);

    static EvalOutT EvalOut;
    bool Feasible;
    POLYQP_TRY(evaluateConstraints(QPTraj, Plan, Constants, Limits, EvalOut, Feasible));
    Real bestUtil = EvalOut.MaxOverallUtil;
    RowSeg bestSegT = Plan.SegT;
    RowSeg prevSegUtil;
    maxDim1(EvalOut.SegUtil, prevSegUtil);

    while (outerIter < MaxOuter) {
        bool didScale = false;
        bool didContract = false;

        for (int leg = 1; leg <= nLegs; ++leg) {
            const int m1 = V1(Plan.LegStartSeg, leg);
            const int m2 = V1(Plan.LegEndSeg, leg);
            const Real s_leg = maxv(EvalOut.SegDilate(R1(m1, m2)));
            if (s_leg > Real(1.0)) {
                const Real s_apply = mmin(Real(2.0), Real(1.02) * s_leg);
                Plan.SegT(R1(m1, m2)) = Plan.SegT(R1(m1, m2)) * s_apply;
                V1(legHasDilated, leg) = true;
                didScale = true;
            } else if (Feasible && V1(legCanContract, leg) && !V1(legHasDilated, leg) && (Real(1.02) * s_leg < Real(0.97))) {
                const Real s_min_leg = maxv((Limits.MinSegTime / Plan.SegT(R1(m1, m2)).array()).matrix());
                const Real s_apply = maxv({Real(0.60), s_min_leg, Real(1.02) * s_leg});
                if (s_apply < Real(1.0) - Real(1e-4)) {
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
        POLYQP_TRY(evaluateConstraints(QPTraj, Plan, Constants, Limits, EvalOut, Feasible));

        if (didContract) {
            bestUtil = EvalOut.MaxOverallUtil;
            if (Feasible) {
                bestSegT = Plan.SegT;
            }
            maxDim1(EvalOut.SegUtil, prevSegUtil);
            stallCount = 0;
            continue;
        }

        // Stall detection during dilation phase
        RowSeg curSegUtil;
        maxDim1(EvalOut.SegUtil, curSegUtil);
        const BoolSeg violMask = curSegUtil.array() > Real(1.0) + Real(1e-3);
        bool anyViolImproved = false;
        for (int q = 1; q <= static_cast<int>(violMask.size()); ++q) {
            if (V1(violMask, q) && V1(curSegUtil, q) < V1(prevSegUtil, q) - Real(1e-4)) anyViolImproved = true;
        }
        const bool madeProgress = (EvalOut.MaxOverallUtil < bestUtil - Real(1e-4)) || anyViolImproved;
        prevSegUtil = curSegUtil;

        if (EvalOut.MaxOverallUtil < bestUtil - Real(1e-4)) {
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
                POLYQP_TRY(evaluateConstraints(QPTraj, Plan, Constants, Limits, EvalOut, Feasible));
                break;
            }
        }
    }

    // 5. Sample uniform output grid and build Report
    POLYQP_TRY(sampleOutputs(QPTraj, Plan, Constants, Limits, InitialGuess));
    POLYQP_TRY(buildReport(Feasible, EvalOut.MaxUtilPerType, outerIter, QPTraj, Report));
    return PolyQPStatus::Ok;
}

static PolyQPStatus solveAllLegs(PlanT &Plan, const MatNC &A0, const MatNC &H0, const LimitsT &Limits, bool doAllocateTimes, QPTrajT &QPTraj) {
    // Solve All Flight Legs with P-Space Fixed-Point Refinement
    const int MAX_SEG = POLYQP_MAX_SEG;
    const int MAX_KNOTS = MAX_SEG + 1;
    const int r = Plan.R;
    const int nc = 2 * r;
    const int M = Plan.SegCount;
    const int K = Plan.KnotCount;
    const Real rho = Limits.Rho;
    const Real g = Limits.g;

    // Initialize P-space boundary conditions: P = r_cg + rho * b3 (init b3 = e3)
    // TRANSPILE-FLAG F03 (fixed by FIX-F03): the knot-1 P-offset assumed b3 = e3 (rho ~ 0.22 m for ASTRA) and the velocity BC was applied to dP/dt, not v_cg + rho*db3/dt; a knot with B3BC/B3DotBC now uses the given attitude.
    static Fix3 FixP;
    sub3(Plan.FixCg, 1, 3, 1, r, 1, K, FixP);
    for (int k = 1; k <= K; ++k) {
        const bool hasB3 = (!Col1(Plan.B3BC, k).array().isNaN()).all(); // FIX-F03: a knot with a given attitude uses P = r_cg + rho * b3
        const bool hasB3Dot = (!Col1(Plan.B3DotBC, k).array().isNaN()).all(); // FIX-F03: and dP/dt = v_cg + rho * db3/dt where the velocity is fixed
        for (int d = 1; d <= 3; ++d) {
            if (!std::isnan(FixP(d, 1, k))) {
                if (hasB3) { // FIX-F03: P = r_cg + rho * b3 with the given attitude
                    FixP(d, 1, k) = Plan.FixCg(d, 1, k) + rho * M1(Plan.B3BC, d, k);
                } else if (d == 3) {
                    FixP(d, 1, k) = Plan.FixCg(d, 1, k) + rho;
                } else {
                    FixP(d, 1, k) = Plan.FixCg(d, 1, k);
                }
            }
            if (hasB3Dot && !std::isnan(FixP(d, 2, k))) { // FIX-F03: dP/dt = v_cg + rho * db3/dt
                FixP(d, 2, k) = Plan.FixCg(d, 2, k) + rho * M1(Plan.B3DotBC, d, k);
            }
        }
    }

    static Coef3 CoefAll;
    CoefAll.setConstant(POLYQP_MAX_NC, MAX_SEG, 3, Real(0));
    static Fix3 DerivAll;
    DerivAll = FixP;
    Real totalCost = 0;

    // Optional one-time Richter time redistribution on multi-segment legs where intermediate knots have at least one free velocity axis
    if (doAllocateTimes) {
        for (int leg = 1; leg <= Plan.LegCount; ++leg) {
            const int m1 = V1(Plan.LegStartSeg, leg);
            const int m2 = V1(Plan.LegEndSeg, leg);
            if (m2 > m1) {
                static Fix3 midVelBins;
                sub3(FixP, 1, 3, 2, 2, m1 + 1, m2, midVelBins);
                Bool3 isnanMidVelBins;
                isnan3(midVelBins, isnanMidVelBins);
                Bool3 anyIsnanMidVelBins;
                any3(isnanMidVelBins, 1, anyIsnanMidVelBins);
                Bool3 allAnyIsnanMidVelBins;
                all3(anyIsnanMidVelBins, allAnyIsnanMidVelBins);
                if (ifTrue(allAnyIsnanMidVelBins)) {
                    static Fix3 FixLeg;
                    sub3(FixP, 1, 3, 1, r, m1, m2 + 1, FixLeg);
                    const RowSeg T_leg = Plan.SegT(R1(m1, m2));
                    RowSeg T_opt;
                    Real CostIgnored;
                    POLYQP_TRY(allocateTimes(FixLeg, T_leg, A0, H0, Limits.MinSegTime, T_opt, CostIgnored));
                    Plan.SegT(R1(m1, m2)) = T_opt;
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
            static Fix3 FixLeg;
            sub3(FixP, 1, 3, 1, r, m1, m2 + 1, FixLeg);
            const RowSeg T_leg = Plan.SegT(R1(m1, m2));

            static Coef3 CoefLeg;
            static Fix3 DerivLeg;
            Real costLeg;
            POLYQP_TRY(solveLeg(FixLeg, T_leg, A0, H0, CoefLeg, DerivLeg, costLeg));
            setSub3(CoefAll, 1, nc, m1, m2, 1, 3, CoefLeg);
            setSub3(DerivAll, 1, 3, 1, r, m1, m2 + 1, DerivLeg);
            totalCost = totalCost + costLeg;
        }

        // Populate Hold and Cruise segment polynomial coefficients
        for (int m = 1; m <= M; ++m) {
            if (V1(Plan.IsHold, m)) {
                for (int d = 1; d <= 3; ++d) {
                    fillSub3(CoefAll, 1, nc, m, m, d, d, Real(0));
                    CoefAll(1, m, d) = FixP(d, 1, m);
                }
            } else if (V1(Plan.IsCruise, m)) {
                for (int d = 1; d <= 3; ++d) {
                    fillSub3(CoefAll, 1, nc, m, m, d, d, Real(0));
                    CoefAll(1, m, d) = FixP(d, 1, m);
                    CoefAll(2, m, d) = FixP(d, 2, m);
                }
            }
        }

        // Update P = r_cg + rho * b3 at pass-through knots using solved acceleration
        if (pass < nPasses) {
            // TRANSPILE-FLAG F03 (fixed by FIX-F03): refinement covers k = 2:(K-1) only; the knot-1 attitude now comes from B3BC instead of the b3 = e3 assumption.
            for (int k = 2; k <= K - 1; ++k) {
                if (!V1(Plan.KnotIsRest, k) && Col1(Plan.B3BC, k).array().isNaN().any()) { // FIX-F03: a given knot attitude is not overwritten by the refinement
                    const Vec3 a_k = getFiber3(DerivAll, 3, k);
                    const Vec3 e_k = a_k + Vec3(0, 0, g);
                    const Real ne = mnorm(e_k);
                    if (ne > Real(1e-6)) {
                        const Vec3 b3_k = e_k / ne;
                        for (int d = 1; d <= 3; ++d) {
                            if (!std::isnan(Plan.FixCg(d, 1, k))) {
                                FixP(d, 1, k) = Plan.FixCg(d, 1, k) + rho * V1(b3_k, d);
                            }
                        }
                    }
                }
            }
        }
    }

    RowKnot KeyTimes = RowKnot::Zero(MAX_KNOTS);
    RowSeg cumsumSegT;
    cumsum(Plan.SegT(R1(1, M)), cumsumSegT);
    KeyTimes(R1(1, M + 1)) << Real(0), cumsumSegT;

    static FixOut3 FixOut;
    FixOut.setConstant(3, 6, MAX_KNOTS, NaN());
    setSub3(FixOut, 1, 3, 1, r, 1, K, DerivAll);

    QPTraj.R = r;
    QPTraj.Count = M;
    QPTraj.T = Plan.SegT;
    QPTraj.KeyTimes = KeyTimes;
    sub3(CoefAll, 1, nc, 1, CoefAll.size(2), 1, CoefAll.size(3), QPTraj.Coef);
    QPTraj.IsHold = Plan.IsHold;
    QPTraj.IsCruise = Plan.IsCruise;
    QPTraj.Fix = FixOut;
    QPTraj.Cost = totalCost;
    QPTraj.T_total = sumv(Plan.SegT(R1(1, M)));
    QPTraj.FlatOffset = rho;
    return PolyQPStatus::Ok;
}

static PolyQPStatus buildReport(bool Feasible, const Vec9 &MaxUtilPerType, int outerIter, const QPTrajT &QPTraj, ReportT &Report) {
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
    Report.SegTimes = QPTraj.T(R1(1, M));
    Report.T_total = QPTraj.T_total;
    return PolyQPStatus::Ok;
}

}
