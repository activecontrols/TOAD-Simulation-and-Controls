// Transpiled from TOAD/Guidance/PolynomialQP/flatnessInvert.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

PolyQPStatus flatnessInvert(const QPTrajT &QPTraj, const RowGrid &t_vec, const Vec3 &HDGRef, const PolyQPConstants &Constants, const LimitsT &Limits,
                            Mat15Grid &X, WrenchT &Wrench, Mat3Grid &AccelCg) {
    // Differential Flatness Inversion at Centre of Percussion
    const int nPts = static_cast<int>(t_vec.size());
    static Mat3Grid P_pos;
    POLYQP_TRY(PolyEval(QPTraj, t_vec, P_pos, 0));
    static Mat3Grid P_vel;
    POLYQP_TRY(PolyEval(QPTraj, t_vec, P_vel, 1));
    static Mat3Grid P_acc;
    POLYQP_TRY(PolyEval(QPTraj, t_vec, P_acc, 2));
    static Mat3Grid P_jerk;
    POLYQP_TRY(PolyEval(QPTraj, t_vec, P_jerk, 3));
    static Mat3Grid P_snap;
    POLYQP_TRY(PolyEval(QPTraj, t_vec, P_snap, 4));

    const Real g = Limits.g;
    const Vec3 G_vec(0, 0, g);
    const Real rho = QPTraj.FlatOffset;

    // Heading reference vector in horizontal plane
    Vec3 x_C(V1(HDGRef, 1), V1(HDGRef, 2), 0);
    if (mnorm(x_C) < Real(1e-9)) {
        x_C = Vec3(1, 0, 0);
    } else {
        x_C = x_C / mnorm(x_C);
    }

    // First pass: Attitude, angular velocity, angular acceleration, and specific force
    static MatX<4, POLYQP_MAX_NODES> q_s;
    q_s.setZero(4, nPts);
    static Mat3Grid b1_s;
    b1_s.setZero(3, nPts);
    static Mat3Grid b2_s;
    b2_s.setZero(3, nPts);
    static Mat3Grid b3_s;
    b3_s.setZero(3, nPts);
    static Mat3Grid b3d_s;
    b3d_s.setZero(3, nPts);
    static Mat3Grid b3dd_s;
    b3dd_s.setZero(3, nPts);
    static Mat3Grid w_s;
    w_s.setZero(3, nPts);
    static Mat3Grid wd_s;
    wd_s.setZero(3, nPts);
    static RowGrid n_s;
    n_s.setZero(nPts);

    for (int k = 1; k <= nPts; ++k) {
        Vec3 e = Col1(P_acc, k) + G_vec;
        const Vec3 ed = Col1(P_jerk, k);
        const Vec3 edd = Col1(P_snap, k);

        Real n = mnorm(e);
        if (n < Real(1e-6)) {
            e = Vec3(0, 0, Real(1e-6));
            n = Real(1e-6);
        }
        const Vec3 u = e / n;
        const Real nd = mdot(u, ed);
        const Vec3 ud = (ed - u * nd) / n;
        const Real ndd = mdot(ud, ed) + mdot(u, edd);
        const Vec3 udd = (edd - 2 * ud * nd - u * ndd) / n;

        // Body triad [b1, b2, b3] with b3 = u and b2 = norm(b3 x x_C)
        Vec3 x_ref = x_C; // FIX-F27: the heading vector actually used for b2, reused for its derivatives
        Vec3 c = mcross(u, x_C);
        Real nc = mnorm(c);
        if (nc < Real(1e-6)) {
            // TRANSPILE-FLAG F27 (fixed by FIX-F27): the x_C_alt fallback changed b2 but cd/cdd still differentiated cross(u, x_C).
            const Vec3 x_C_alt(0, 1, 0);
            c = mcross(u, x_C_alt);
            nc = mnorm(c);
            x_ref = x_C_alt; // FIX-F27: the fallback heading also drives b2d and b2dd
        }
        const Vec3 b2 = c / nc;
        const Vec3 b1 = mcross(b2, u);

        // First and second time derivatives of b2 and b1
        const Vec3 cd = mcross(ud, x_ref); // FIX-F27: was cross(ud, x_C) even when the fallback heading was used
        const Real ncd = mdot(b2, cd);
        const Vec3 b2d = (cd - b2 * ncd) / nc;
        const Vec3 cdd = mcross(udd, x_ref); // FIX-F27: was cross(udd, x_C) even when the fallback heading was used
        const Real ncdd = mdot(b2d, cd) + mdot(b2, cdd);
        const Vec3 b2dd = (cdd - 2 * b2d * ncd - b2 * ncdd) / nc;

        const Vec3 b1d = mcross(b2d, u) + mcross(b2, ud);
        const Vec3 b1dd = mcross(b2dd, u) + 2 * mcross(b2d, ud) + mcross(b2, udd);

        // Body angular velocity w_B and angular acceleration wd_B
        const Real wx = -mdot(ud, b2);
        const Real wy = mdot(ud, b1);
        const Real wz = mdot(b1d, b2);

        const Real wxd = -mdot(udd, b2) - mdot(ud, b2d);
        const Real wyd = mdot(udd, b1) + mdot(ud, b1d);
        const Real wzd = mdot(b1dd, b2) + mdot(b1d, b2d);

        // Quaternion from DCM C_IB = [b1, b2, b3] (DCM_Quat_Conversion expects C_IB)
        Mat3 C_IB;
        C_IB << b1, b2, u;
        Vec4 qk = DCM_Quat_Conversion(C_IB);
        if (k > 1 && mdot(Col1(q_s, k - 1), qk) < 0) {
            qk = -qk;
        }

        Col1(q_s, k) = qk;
        Col1(b1_s, k) = b1;
        Col1(b2_s, k) = b2;
        Col1(b3_s, k) = u;
        Col1(b3d_s, k) = ud;
        Col1(b3dd_s, k) = udd;
        Col1(w_s, k) = Vec3(wx, wy, wz);
        Col1(wd_s, k) = Vec3(wxd, wyd, wzd);
        V1(n_s, k) = n;
    }

    // CG kinematics from P-space: r_cg = P - rho * b3
    static Mat3Grid r_s;
    r_s = P_pos - rho * b3_s;
    static Mat3Grid v_s;
    v_s = P_vel - rho * b3d_s;
    AccelCg = P_acc - rho * b3dd_s;

    // Second pass: Mass depletion and wrench synthesis
    static RowGrid ml_s;
    ml_s.setZero(nPts);
    static RowGrid mi_s;
    mi_s.setZero(nPts);
    if (Limits.Vehicle == 1) {
        // TRANSPILE-FLAG F29 (fixed by FIX-F29): propellant integration started from the full-tank masses; Limits.OxMass/FuMass are now the estimate when Opts.InitOxMass/InitFuMass are given (buildVehicleLimits).
        V1(ml_s, 1) = Limits.OxMass;
        V1(mi_s, 1) = Limits.FuMass;
    }

    static Mat3Grid F_body;
    F_body.setZero(3, nPts);
    static Mat3Grid Tau_body;
    Tau_body.setZero(3, nPts);
    static RowGrid CGz_vec;
    CGz_vec.setZero(nPts);

    for (int k = 1; k <= nPts; ++k) {
        const Real m_k = Limits.m_dry + V1(ml_s, k) + V1(mi_s, k);
        Mat3 J_tot;
        Real CGz;
        if (Limits.HasVehicle) {
            ComputeJtot(V1(ml_s, k), V1(mi_s, k), Constants, J_tot, CGz);
        } else {
            J_tot = Limits.J_init;
            CGz = Limits.CGz_init;
        }
        V1(CGz_vec, k) = CGz;

        // Rotational Euler dynamics: M_B = J_tot * wd + w x (J_tot * w)
        const Vec3 wk = Col1(w_s, k);
        const Vec3 wdk = Col1(wd_s, k);
        const Vec3 tauk = J_tot * wdk + mcross(wk, J_tot * wk);
        Col1(Tau_body, k) = tauk;

        // Exact body force required to track CG acceleration: F_B = m * C_BI * (a_cg + G)
        Mat3 C_BI;
        C_BI << Col1(b1_s, k).transpose(), Col1(b2_s, k).transpose(), Col1(b3_s, k).transpose();
        const Vec3 C_BI_aG = C_BI * (Col1(AccelCg, k) + G_vec);
        const Vec3 Fk = m_k * C_BI_aG;
        Col1(F_body, k) = Fk;

        // Forward integrate propellant depletion for Vehicle 1
        if (k < nPts && Limits.Vehicle == 1) {
            const Real dt_k = V1(t_vec, k + 1) - V1(t_vec, k);
            const Real T_mag = mnorm(Fk);
            const Real OF = Limits.OF;
            const Real mdot_tot = (T_mag / Limits.MaxThrust) * Limits.MaxMdot;
            V1(ml_s, k + 1) = mmax(Real(0), V1(ml_s, k) - dt_k * (OF / (1 + OF)) * mdot_tot);
            V1(mi_s, k + 1) = mmax(Real(0), V1(mi_s, k) - dt_k * (1 / (1 + OF)) * mdot_tot);
        }
    }

    X.resize(15, nPts);
    X << q_s, r_s, v_s, w_s, ml_s, mi_s;
    Wrench.F_body = F_body;
    Wrench.Tau_body = Tau_body;
    Wrench.CGz = CGz_vec;
    return PolyQPStatus::Ok;
}

}
