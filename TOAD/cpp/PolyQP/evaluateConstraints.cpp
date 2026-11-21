// Transpiled from TOAD/Guidance/PolynomialQP/evaluateConstraints.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

static Real solveThrustStretch(const Vec3 &a, const Vec3 &G_vec, Real e_allow, Real u_thi);
static Real solveThrustFloorStretch(const Vec3 &a, const Vec3 &G_vec, Real e_floor, Real u_tlo);

PolyQPStatus evaluateConstraints(const QPTrajT &QPTraj, const PlanT &Plan, const PolyQPConstants &Constants, const LimitsT &Limits, EvalOutT &EvalOut,
                                 bool &Feasible) {
    // Dense-Grid Constraint Utilisation & Per-Segment Dilation Ratios
    const int M = QPTraj.Count;
    const Real T_total = QPTraj.T_total;
    const Real dt_eval = Real(0.5) * Limits.DtSample;
    // TRANSPILE-FLAG F11: the 3999 cap silently truncates the evaluation grid for long trajectories. Transpiled as-is per MATLAB.
    const int nGrid = static_cast<int>(mmin(Real(3999), mmax(Real(41), std::ceil(T_total / dt_eval) + 1)));
    POLYQP_CAPACITY(nGrid <= POLYQP_MAX_GRID);
    static RowGrid t_grid;
    linspace(0, T_total, nGrid, t_grid);
    const Real dt_step = V1(t_grid, 2) - V1(t_grid, 1);

    // Invert flat trajectory and allocate controls across grid
    static Mat15Grid X;
    static WrenchT Wrench;
    static Mat3Grid AccelCg;
    POLYQP_TRY(flatnessInvert(QPTraj, t_grid, Plan.HDGRef, Constants, Limits, X, Wrench, AccelCg));
    static Mat4Grid U;
    POLYQP_TRY(allocateControls(Wrench, Limits, U));

    // Finite-difference actuator rates on dense grid
    static Mat4Grid U_rate;
    U_rate.setZero(4, nGrid);
    if (nGrid >= 2) {
        static Mat4Grid diffU;
        diff_dim2(U, diffU);
        U_rate(Eigen::all, R1(1, nGrid - 1)) = diffU / dt_step;
        Col1(U_rate, nGrid) = Col1(U_rate, nGrid - 1);
    }

    // Map grid samples to segments
    const RowKnot kt = QPTraj.KeyTimes(R1(1, M + 1));
    static IntGrid seg_idx;
    seg_idx.setConstant(nGrid, 1);
    for (int k = 1; k <= nGrid; ++k) {
        const Real tk = V1(t_grid, k);
        int m_found = 1;
        for (int m = 1; m <= M; ++m) {
            if (tk > V1(kt, m)) {
                m_found = m;
            } else if (tk == V1(kt, m) && !V1(QPTraj.IsHold, m) && !V1(Plan.IsCruise, m)) {
                m_found = m;
            }
        }
        V1(seg_idx, k) = m_found;
    }

    // FIX-F01: first grid sample of each segment, whose CG velocity is the plan's start velocity for D3
    IntSeg kFirst = IntSeg::Zero(M);
    for (int k = nGrid; k >= 1; --k) {
        V1(kFirst, V1(seg_idx, k)) = k;
    }

    // Zero finite-difference actuator rates at hold/cruise segment boundaries
    for (int k = 1; k <= nGrid - 1; ++k) {
        const int m1 = V1(seg_idx, k);
        const int m2 = V1(seg_idx, k + 1);
        if (V1(QPTraj.IsHold, m1) || V1(QPTraj.IsHold, m2) || V1(Plan.IsCruise, m1) || V1(Plan.IsCruise, m2)) {
            Col1(U_rate, k).setZero();
        }
    }
    Col1(U_rate, nGrid) = Col1(U_rate, std::max(1, nGrid - 1));

    // Utilisation arrays per segment:
    // Rows: 1=Vel, 2=Accel, 3=ThrustHi, 4=ThrustLo, 5=Gimbal, 6=GimbalRate, 7=ThrustRate, 8=Roll, 9=Tilt
    Mat9Seg SegUtil = Mat9Seg::Zero(9, M);
    RowSeg SegDilate = RowSeg::Zero(M);

    const Real g = Limits.g;
    const Vec3 G_vec(0, 0, g);

    for (int k = 1; k <= nGrid; ++k) {
        const int m = V1(seg_idx, k);
        if (V1(QPTraj.IsHold, m)) {
            continue;
        }

        const Vec3 vk = mabs(Col1(X, k)(R1(8, 10)));
        const Vec3 ak = Col1(AccelCg, k);
        const Vec4 qk = Col1(X, k)(R1(1, 4));
        const Vec4 uk = Col1(U, k);
        const Vec4 urk = mabs(Col1(U_rate, k));
        const Real mk = Limits.m_dry + M1(X, 14, k) + M1(X, 15, k);

        // 1. Velocity utilisation (k = 1)
        // FIX-F01: a pinned start-knot velocity raises that axis's cap to the larger of |v_BC| and the plan's CG speed at the segment start (D3)
        Vec3 vCap = Col1(Plan.SegMaxVel, m);
        const Vec3 vBC = mabs(getFiber3(Plan.FixCg, 2, m));
        Vec3 vStart = Vec3::Zero();
        if (V1(kFirst, m) > 0) {
            vStart = mabs(Vec3(Col1(X, V1(kFirst, m))(R1(8, 10))));
        }
        for (int d = 1; d <= 3; ++d) {
            if (!std::isnan(V1(vBC, d))) {
                V1(vCap, d) = maxv({V1(vCap, d), V1(vBC, d), V1(vStart, d)});
            }
        }
        const Real u_vel = maxv((vk.array() / vCap.array()).matrix()); // FIX-F01: cap raised per D3

        // 2. Acceleration utilisation (k = 2)
        const Real u_acc = maxv((mabs(ak).array() / Col1(Plan.SegMaxAccel, m).array()).matrix());

        // 3. Thrust upper & lower utilisation
        const Real T_k = V1(uk, 3);
        const Real u_thi = T_k / Limits.T_hi;
        const Real u_tlo = Limits.T_lo / mmax(Real(1e-3), T_k);

        // 4. Gimbal angle (k = 4)
        const Real u_gim = maxv(mabs(uk(R1(1, 2)))) / Limits.GimbalLimit;

        // 5. Gimbal rate (k = 5)
        Real u_grate = maxv(urk(R1(1, 2))) / Limits.GimbalRateLim;

        // 6. Thrust rate (k = 3)
        Real u_trate = V1(urk, 3) / Limits.ThrustRateLim;

        // 7. Roll torque & rate (k = 4)
        const Real u_roll = mmax(std::fabs(V1(uk, 4)) / Limits.RollLimit, V1(urk, 4) / Limits.RollRateLim);

        // 8. Attitude tilt angle from vertical (k = 2)
        const Real b3_z = 1 - 2 * (std::pow(V1(qk, 2), Real(2)) + std::pow(V1(qk, 3), Real(2)));
        const Real tilt_k = std::acos(mmax(Real(-1), mmin(Real(1), b3_z)));
        const Real u_tilt = tilt_k / Limits.MaxTilt;

        // Skip recording finite-difference actuator rates on exact 1D cruise segments
        if (V1(Plan.IsCruise, m)) {
            u_grate = 0;
            u_trate = 0;
        }

        Vec9 u_col(9);
        u_col << u_vel, u_acc, u_thi, u_tlo, u_gim, u_grate, u_trate, u_roll, u_tilt;
        Col1(SegUtil, m) = mmax(Col1(SegUtil, m), u_col);

        // Compute bidirectional segment time-scale factor sigma_m = max_c (u_c)^(1/k_c)
        if (!V1(Plan.IsCruise, m)) {
            const Real s_thi = solveThrustStretch(ak, G_vec, Limits.T_hi / mk, u_thi);
            const Real s_tlo = solveThrustFloorStretch(ak, G_vec, Limits.T_lo / mk, u_tlo);

            // TRANSPILE-FLAG F12: dilation exponents differ from paper Sec. 3.5.2 (T^-3); the code is followed. Transpiled as-is per MATLAB.
            const Real s_req = maxv({u_vel, std::pow(u_acc, Real(1) / Real(2)), s_thi, s_tlo, std::pow(u_gim, Real(1) / Real(4)),
                                     std::pow(u_grate, Real(1) / Real(5)), std::pow(u_trate, Real(1) / Real(3)), std::pow(u_roll, Real(1) / Real(4)),
                                     std::pow(u_tilt, Real(1) / Real(2))});

            V1(SegDilate, m) = mmax(V1(SegDilate, m), s_req);
        }
    }

    for (int m = 1; m <= M; ++m) {
        if (V1(QPTraj.IsHold, m) || V1(Plan.IsCruise, m) || V1(SegDilate, m) <= 0) {
            V1(SegDilate, m) = Real(1.0);
        }
    }

    Vec9 maxUtilPerType;
    maxDim2(SegUtil, maxUtilPerType);
    const Real maxOverallUtil = maxv(maxUtilPerType);
    Feasible = (maxOverallUtil <= Real(1.0) + Real(1e-3));

    EvalOut.SegUtil = SegUtil;
    EvalOut.MaxUtilPerType = maxUtilPerType;
    EvalOut.MaxOverallUtil = maxOverallUtil;
    EvalOut.SegDilate = SegDilate;
    return PolyQPStatus::Ok;
}

static Real solveThrustStretch(const Vec3 &a, const Vec3 &G_vec, Real e_allow, Real u_thi) {
    // Quadratic root for time scale s where specific force reaches e_allow
    const Real g = V1(G_vec, 3);
    const Real a2 = mdot(a, a);
    const Real aG = mdot(a, G_vec);
    const Real c = std::pow(g, Real(2)) - std::pow(e_allow, Real(2));
    Real s;
    if (c >= 0) {
        s = std::sqrt(u_thi);
        return s;
    }
    if (a2 < Real(1e-8)) {
        s = Real(0.0);
        return s;
    }
    const Real disc = mmax(Real(0), std::pow(aG, Real(2)) - a2 * c);
    const Real w = (-aG + std::sqrt(disc)) / a2;
    if (w > 0) {
        s = 1 / std::sqrt(w);
    } else {
        s = std::sqrt(u_thi);
    }
    return s;
}

static Real solveThrustFloorStretch(const Vec3 &a, const Vec3 &G_vec, Real e_floor, Real u_tlo) {
    // Time scale s where downward acceleration causes thrust to reach T_lo
    const Real g = V1(G_vec, 3);
    const Real a2 = mdot(a, a);
    const Real aG = mdot(a, G_vec);
    const Real c = std::pow(g, Real(2)) - std::pow(e_floor, Real(2));
    Real s;
    if (c <= 0) {
        // Hover thrust m*g is already below T_lo: time dilation cannot raise thrust
        s = Real(1.0);
        return s;
    }
    if (a2 < Real(1e-8) || aG >= 0) {
        s = Real(0.0);
        return s;
    }
    const Real disc = std::pow(aG, Real(2)) - a2 * c;
    if (disc < 0) {
        if (u_tlo > Real(1.0)) {
            s = Real(1.0);
        } else {
            s = Real(0.0);
        }
        return s;
    }
    const Real w = (-aG - std::sqrt(disc)) / a2;
    if (w > 0) {
        s = 1 / std::sqrt(w);
    } else {
        s = Real(0.0);
    }
    return s;
}

}
