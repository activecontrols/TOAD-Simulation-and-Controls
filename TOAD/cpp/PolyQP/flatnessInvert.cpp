// Transpiled from TOAD/Guidance/PolynomialQP/flatnessInvert.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

PolyQPStatus flatnessInvertPoint(const QPTrajT &QPTraj, double t, const Vector3d &HDGRef, const LimitsT &Limits, bool hasPrevQ, Vector4d &qPrev, FlatPointT &Out) {
    // Differential Flatness Inversion at Centre of Percussion
    // STREAM: flatnessInvert.m's loops for one time sample. Without the vehicle 1 mass integration, every sample depends only on
    // STREAM: its own time, except the quaternion sign continuity, which needs the previous sample's quaternion (qPrev).
    const Vector3d P_pos = PolyEvalPoint(QPTraj, t, 0);
    const Vector3d P_vel = PolyEvalPoint(QPTraj, t, 1);
    const Vector3d P_acc = PolyEvalPoint(QPTraj, t, 2);
    const Vector3d P_jerk = PolyEvalPoint(QPTraj, t, 3);
    const Vector3d P_snap = PolyEvalPoint(QPTraj, t, 4);

    const double g = Limits.g;
    const Vector3d G_vec(0, 0, g);
    const double rho = QPTraj.FlatOffset;

    // Heading reference vector in horizontal plane
    Vector3d x_C(V1(HDGRef, 1), V1(HDGRef, 2), 0);
    if (mnorm(x_C) < double(1e-9)) {
        x_C = Vector3d(1, 0, 0);
    } else {
        x_C = x_C / mnorm(x_C);
    }

    // First pass: Attitude, angular velocity, angular acceleration, and specific force
    Vector3d e = P_acc + G_vec;
    const Vector3d ed = P_jerk;
    const Vector3d edd = P_snap;

    double n = mnorm(e);
    if (n < double(1e-6)) {
        e = Vector3d(0, 0, double(1e-6));
        n = double(1e-6);
    }
    const Vector3d u = e / n;
    const double nd = mdot(u, ed);
    const Vector3d ud = (ed - u * nd) / n;
    const double ndd = mdot(ud, ed) + mdot(u, edd);
    const Vector3d udd = (edd - 2 * ud * nd - u * ndd) / n;

    // Body triad [b1, b2, b3] with b3 = u and b2 = norm(b3 x x_C)
    Vector3d x_ref = x_C; // FIX-F27: the heading vector actually used for b2, reused for its derivatives
    Vector3d c = mcross(u, x_C);
    double nc = mnorm(c);
    if (nc < double(1e-6)) {
        // TRANSPILE-FLAG F27 (fixed by FIX-F27): the x_C_alt fallback changed b2 but cd/cdd still differentiated cross(u, x_C).
        const Vector3d x_C_alt(0, 1, 0);
        c = mcross(u, x_C_alt);
        nc = mnorm(c);
        x_ref = x_C_alt; // FIX-F27: the fallback heading also drives b2d and b2dd
    }
    const Vector3d b2 = c / nc;
    const Vector3d b1 = mcross(b2, u);

    // First and second time derivatives of b2 and b1
    const Vector3d cd = mcross(ud, x_ref); // FIX-F27: was cross(ud, x_C) even when the fallback heading was used
    const double ncd = mdot(b2, cd);
    const Vector3d b2d = (cd - b2 * ncd) / nc;
    const Vector3d cdd = mcross(udd, x_ref); // FIX-F27: was cross(udd, x_C) even when the fallback heading was used
    const double ncdd = mdot(b2d, cd) + mdot(b2, cdd);
    const Vector3d b2dd = (cdd - 2 * b2d * ncd - b2 * ncdd) / nc;

    const Vector3d b1d = mcross(b2d, u) + mcross(b2, ud);
    const Vector3d b1dd = mcross(b2dd, u) + 2 * mcross(b2d, ud) + mcross(b2, udd);

    // Body angular velocity w_B and angular acceleration wd_B
    const double wx = -mdot(ud, b2);
    const double wy = mdot(ud, b1);
    const double wz = mdot(b1d, b2);

    const double wxd = -mdot(udd, b2) - mdot(ud, b2d);
    const double wyd = mdot(udd, b1) + mdot(ud, b1d);
    const double wzd = mdot(b1dd, b2) + mdot(b1d, b2d);

    // Quaternion from DCM C_IB = [b1, b2, b3] (DCM_Quat_Conversion expects C_IB)
    Matrix3_3d C_IB;
    C_IB << b1, b2, u;
    Vector4d qk = DCM_Quat_Conversion(C_IB);
    if (hasPrevQ && mdot(qPrev, qk) < 0) {
        qk = -qk;
    }
    qPrev = qk;

    // CG kinematics from P-space: r_cg = P - rho * b3
    const Vector3d r_s = P_pos - rho * u;
    const Vector3d v_s = P_vel - rho * ud;
    const Vector3d AccelCg = P_acc - rho * udd;

    // Second pass: Mass depletion and wrench synthesis
    // ASTRA-ONLY: ASTRA carries no propellant, so ml_s and mi_s stay zero (X(14:15) keep their slots); the vehicle 1 propellant-depletion integration was removed.
    const double ml_s = 0;
    const double mi_s = 0;

    const double m_k = Limits.m_dry + ml_s + mi_s;
    Matrix3_3d J_tot;
    double CGz;
    if (Limits.HasVehicle) {
        ComputeJtot(ml_s, mi_s, J_tot, CGz);
    } else {
        J_tot = Limits.J_init;
        CGz = Limits.CGz_init;
    }

    // Rotational Euler dynamics: M_B = J_tot * wd + w x (J_tot * w)
    const Vector3d wk(wx, wy, wz);
    const Vector3d wdk(wxd, wyd, wzd);
    const Vector3d tauk = J_tot * wdk + mcross(wk, J_tot * wk);

    // Exact body force required to track CG acceleration: F_B = m * C_BI * (a_cg + G)
    Matrix3_3d C_BI;
    C_BI << b1.transpose(), b2.transpose(), u.transpose();
    const Vector3d C_BI_aG = C_BI * (AccelCg + G_vec);
    const Vector3d Fk = m_k * C_BI_aG;

    Out.X << qk, r_s, v_s, wk, ml_s, mi_s;
    Out.F_body = Fk;
    Out.Tau_body = tauk;
    Out.CGz = CGz;
    Out.AccelCg = AccelCg;
    return PolyQPStatus::Ok;
}

PolyQPStatus flatnessInvert(const QPTrajT &QPTraj, const Matrix1_Nd &t_vec, int nPts, const Vector3d &HDGRef, const LimitsT &Limits, Matrix15_Nd &X, WrenchT &Wrench,
                            Matrix3_Nd &AccelCg) {
    // STREAM: the grid form (kept for verification): flatnessInvertPoint at each of the nPts times, chaining the quaternion.
    POLYQP_CAPACITY(nPts <= POLYQP_MAX_NODES);
    Vector4d qPrev = Vector4d::Zero();
    FlatPointT pt;
    for (int k = 1; k <= nPts; ++k) {
        POLYQP_TRY(flatnessInvertPoint(QPTraj, V1(t_vec, k), HDGRef, Limits, k > 1, qPrev, pt));
        Col1(X, k) = pt.X;
        Col1(Wrench.F_body, k) = pt.F_body;
        Col1(Wrench.Tau_body, k) = pt.Tau_body;
        V1(Wrench.CGz, k) = pt.CGz;
        Col1(AccelCg, k) = pt.AccelCg;
    }
    return PolyQPStatus::Ok;
}

}
