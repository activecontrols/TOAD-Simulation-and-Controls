// Transpiled from TOAD/Guidance/PolynomialQP/evaluateConstraints.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

static double solveThrustStretch(const Vector3d &a, const Vector3d &G_vec, double e_allow, double u_thi);
static double solveThrustFloorStretch(const Vector3d &a, const Vector3d &G_vec, double e_floor, double u_tlo);

PolyQPStatus evaluateConstraints(const QPTrajT &QPTraj, const PlanT &Plan, const LimitsT &Limits, EvalOutT &EvalOut, bool &Feasible) {
    // Dense-Grid Constraint Utilisation & Per-Segment Dilation Ratios
    // STREAM: the grid is not stored. Each sample is evaluated once (flatnessInvertPoint, allocateControlsPoint); a sample is scored
    // STREAM: one step later, when the next sample's controls give its finite-difference rate. SegUtil and SegDilate accumulate as
    // STREAM: running maxima in grid order, so the numbers are those of the stored-grid MATLAB.
    const int M = QPTraj.Count;
    const double T_total = QPTraj.T_total;
    const double dt_eval = double(0.5) * Limits.DtSample;
    // TRANSPILE-FLAG F11: the 3999 cap silently truncates the evaluation grid for long trajectories. Transpiled as-is per MATLAB.
    const int nGrid = static_cast<int>(std::min(double(3999), std::max(double(41), std::ceil(T_total / dt_eval) + 1)));
    // STREAM: t_grid(k) = linspaceAt(0, T_total, nGrid, k - 1), the element MATLAB's linspace(0, T_total, nGrid) holds
    const double dt_step = linspaceAt(0, T_total, nGrid, 1) - linspaceAt(0, T_total, nGrid, 0);

    // Invert flat trajectory and allocate controls across grid
    // STREAM: per sample; qPrev carries the quaternion sign continuity of flatnessInvert along the grid
    Vector4d qPrev = Vector4d::Zero();
    FlatPointT pt;
    FlatPointT ptPrev;

    // Finite-difference actuator rates on dense grid
    // STREAM: U_rate(:, k) = (U(:, k + 1) - U(:, k)) / dt_step needs the next sample's U; the last sample reuses the previous rate
    Vector4d U_k = Vector4d::Zero();
    Vector4d U_prev = Vector4d::Zero();
    Vector4d U_rate_last = Vector4d::Zero();

    // Map grid samples to segments
    const auto segOf = [&](double tk) {
        int m_found = 1;
        for (int m = 1; m <= M; ++m) {
            if (tk > V1(QPTraj.KeyTimes, m)) {
                m_found = m;
            } else if (tk == V1(QPTraj.KeyTimes, m) && !V1(QPTraj.IsHold, m) && !V1(Plan.IsCruise, m)) {
                m_found = m;
            }
        }
        return m_found;
    };
    int seg_k = 1;
    int seg_prev = 1;

    // FIX-F01: first grid sample of each segment, whose CG velocity is the plan's start velocity for D3
    // STREAM: seg_idx is non-decreasing along the grid, so the first sample of segment m is the first one evaluated in m
    Array1_Si kFirst = Array1_Si::Zero();
    Matrix3_Sd vStartSeg = Matrix3_Sd::Zero(); // STREAM: |v_cg| at kFirst(m), kept instead of the stored X(8:10, kFirst(m))

    // Zero finite-difference actuator rates at hold/cruise segment boundaries
    const auto rateZeroed = [&](int m1, int m2) { return V1(QPTraj.IsHold, m1) || V1(QPTraj.IsHold, m2) || V1(Plan.IsCruise, m1) || V1(Plan.IsCruise, m2); };

    // Utilisation arrays per segment:
    // Rows: 1=Vel, 2=Accel, 3=ThrustHi, 4=ThrustLo, 5=Gimbal, 6=GimbalRate, 7=ThrustRate, 8=Roll, 9=Tilt
    Matrix9_Sd SegUtil = Matrix9_Sd::Zero();
    Matrix1_Sd SegDilate = Matrix1_Sd::Zero();

    const double g = Limits.g;
    const Vector3d G_vec(0, 0, g);

    // STREAM: scores one sample (the body of the MATLAB loop over k), given its finite-difference rate
    const auto scoreSample = [&](const FlatPointT &P, const Vector4d &U_s, int m, const Vector4d &U_rate_s) {
        if (V1(QPTraj.IsHold, m)) {
            return;
        }

        const Vector3d vk = mabs(P.X(R1(8, 10)));
        const Vector3d ak = P.AccelCg;
        const Vector4d qk = P.X(R1(1, 4));
        const Vector4d uk = U_s;
        const Vector4d urk = mabs(U_rate_s);
        const double mk = Limits.m_dry + V1(P.X, 14) + V1(P.X, 15);

        // 1. Velocity utilisation (k = 1)
        // FIX-F01: a pinned start-knot velocity raises that axis's cap to the larger of |v_BC| and the plan's CG speed at the segment start (D3)
        Vector3d vCap = Col1(Plan.SegMaxVel, m);
        const Vector3d vBC = mabs(getFiber3(Plan.FixCg, 2, m));
        const Vector3b vBCFixed = getFiber3(Plan.FixCgFixed, 2, m); // NO-NAN: ~isnan(vBC)
        Vector3d vStart = Vector3d::Zero();
        if (V1(kFirst, m) > 0) {
            vStart = Col1(vStartSeg, m);
        }
        for (int d = 1; d <= 3; ++d) {
            if (V1(vBCFixed, d)) {
                V1(vCap, d) = maxv({V1(vCap, d), V1(vBC, d), V1(vStart, d)});
            }
        }
        const double u_vel = maxv((vk.array() / vCap.array()).matrix()); // FIX-F01: cap raised per D3

        // 2. Acceleration utilisation (k = 2)
        const double u_acc = maxv((mabs(ak).array() / Col1(Plan.SegMaxAccel, m).array()).matrix());

        // 3. Thrust upper & lower utilisation
        const double T_k = V1(uk, 3);
        const double u_thi = T_k / Limits.T_hi;
        const double u_tlo = Limits.T_lo / std::max(double(1e-3), T_k);

        // 4. Gimbal angle (k = 4)
        const double u_gim = maxv(mabs(uk(R1(1, 2)))) / Limits.GimbalLimit;

        // 5. Gimbal rate (k = 5)
        double u_grate = maxv(urk(R1(1, 2))) / Limits.GimbalRateLim;

        // 6. Thrust rate (k = 3)
        double u_trate = V1(urk, 3) / Limits.ThrustRateLim;

        // 7. Roll torque & rate (k = 4)
        const double u_roll = std::max(std::fabs(V1(uk, 4)) / Limits.RollLimit, V1(urk, 4) / Limits.RollRateLim);

        // 8. Attitude tilt angle from vertical (k = 2)
        const double b3_z = 1 - 2 * (std::pow(V1(qk, 2), double(2)) + std::pow(V1(qk, 3), double(2)));
        const double tilt_k = std::acos(std::max(double(-1), std::min(double(1), b3_z)));
        const double u_tilt = tilt_k / Limits.MaxTilt;

        // Skip recording finite-difference actuator rates on exact 1D cruise segments
        if (V1(Plan.IsCruise, m)) {
            u_grate = 0;
            u_trate = 0;
        }

        Vector9d u_col;
        u_col << u_vel, u_acc, u_thi, u_tlo, u_gim, u_grate, u_trate, u_roll, u_tilt;
        Col1(SegUtil, m) = Col1(SegUtil, m).cwiseMax(u_col);

        // Compute bidirectional segment time-scale factor sigma_m = max_c (u_c)^(1/k_c)
        if (!V1(Plan.IsCruise, m)) {
            const double s_thi = solveThrustStretch(ak, G_vec, Limits.T_hi / mk, u_thi);
            const double s_tlo = solveThrustFloorStretch(ak, G_vec, Limits.T_lo / mk, u_tlo);

            // TRANSPILE-FLAG F12: dilation exponents differ from paper Sec. 3.5.2 (T^-3); the code is followed. Transpiled as-is per MATLAB.
            const double s_req = maxv({u_vel, std::pow(u_acc, double(1) / double(2)), s_thi, s_tlo, std::pow(u_gim, double(1) / double(4)),
                                       std::pow(u_grate, double(1) / double(5)), std::pow(u_trate, double(1) / double(3)), std::pow(u_roll, double(1) / double(4)),
                                       std::pow(u_tilt, double(1) / double(2))});

            V1(SegDilate, m) = std::max(V1(SegDilate, m), s_req);
        }
    };

    // STREAM: evaluate sample k, then score sample k - 1 with U_rate(:, k - 1), zeroed at hold/cruise boundaries
    for (int k = 1; k <= nGrid; ++k) {
        const double tk = linspaceAt(0, T_total, nGrid, k - 1);
        POLYQP_TRY(flatnessInvertPoint(QPTraj, tk, Plan.HDGRef, Limits, k > 1, qPrev, pt));
        U_k = allocateControlsPoint(pt.F_body, pt.Tau_body, pt.CGz, Limits);
        seg_k = segOf(tk);
        if (V1(kFirst, seg_k) == 0) {
            V1(kFirst, seg_k) = k;
            Col1(vStartSeg, seg_k) = mabs(pt.X(R1(8, 10)));
        }
        if (k >= 2) {
            Vector4d U_rate_prev = (U_k - U_prev) / dt_step;
            if (rateZeroed(seg_prev, seg_k)) {
                U_rate_prev.setZero();
            }
            U_rate_last = U_rate_prev;
            scoreSample(ptPrev, U_prev, seg_prev, U_rate_prev);
        }
        ptPrev = pt;
        U_prev = U_k;
        seg_prev = seg_k;
    }
    // STREAM: U_rate(:, nGrid) = U_rate(:, nGrid - 1), after the boundary zeroing
    scoreSample(ptPrev, U_prev, seg_prev, U_rate_last);

    for (int m = 1; m <= M; ++m) {
        if (V1(QPTraj.IsHold, m) || V1(Plan.IsCruise, m) || V1(SegDilate, m) <= 0) {
            V1(SegDilate, m) = double(1.0);
        }
    }

    Vector9d maxUtilPerType;
    maxDim2(SegUtil.leftCols(M), maxUtilPerType);
    const double maxOverallUtil = maxv(maxUtilPerType);
    Feasible = (maxOverallUtil <= double(1.0) + double(1e-3));

    EvalOut.SegUtil = SegUtil;
    EvalOut.MaxUtilPerType = maxUtilPerType;
    EvalOut.MaxOverallUtil = maxOverallUtil;
    EvalOut.SegDilate = SegDilate;
    return PolyQPStatus::Ok;
}

static double solveThrustStretch(const Vector3d &a, const Vector3d &G_vec, double e_allow, double u_thi) {
    // Quadratic root for time scale s where specific force reaches e_allow
    const double g = V1(G_vec, 3);
    const double a2 = mdot(a, a);
    const double aG = mdot(a, G_vec);
    const double c = std::pow(g, double(2)) - std::pow(e_allow, double(2));
    double s;
    if (c >= 0) {
        s = std::sqrt(u_thi);
        return s;
    }
    if (a2 < double(1e-8)) {
        s = double(0.0);
        return s;
    }
    const double disc = std::max(double(0), std::pow(aG, double(2)) - a2 * c);
    const double w = (-aG + std::sqrt(disc)) / a2;
    if (w > 0) {
        s = 1 / std::sqrt(w);
    } else {
        s = std::sqrt(u_thi);
    }
    return s;
}

static double solveThrustFloorStretch(const Vector3d &a, const Vector3d &G_vec, double e_floor, double u_tlo) {
    // Time scale s where downward acceleration causes thrust to reach T_lo
    const double g = V1(G_vec, 3);
    const double a2 = mdot(a, a);
    const double aG = mdot(a, G_vec);
    const double c = std::pow(g, double(2)) - std::pow(e_floor, double(2));
    double s;
    if (c <= 0) {
        // Hover thrust m*g is already below T_lo: time dilation cannot raise thrust
        s = double(1.0);
        return s;
    }
    if (a2 < double(1e-8) || aG >= 0) {
        s = double(0.0);
        return s;
    }
    const double disc = std::pow(aG, double(2)) - a2 * c;
    if (disc < 0) {
        if (u_tlo > double(1.0)) {
            s = double(1.0);
        } else {
            s = double(0.0);
        }
        return s;
    }
    const double w = (-aG - std::sqrt(disc)) / a2;
    if (w > 0) {
        s = 1 / std::sqrt(w);
    } else {
        s = double(0.0);
    }
    return s;
}

}
