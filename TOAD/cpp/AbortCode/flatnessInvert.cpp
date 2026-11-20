// Differential Flatness Inversion at CoP
#include "AbortCode.h"
#include <cmath>

void flatnessInvertPoint(const Vector3 &P_pos, const Vector3 &P_vel, const Vector3 &P_acc, const Vector3 &P_jerk, const Vector3 &P_snap, float rho, const constantsASTRA_t &constants,
                         const Vector4 *q_prev, Vector15 &X_out, AbortWrenchPoint &Wrench_out, Vector3 &AccelCg_out) {

    Vector3 G_vec(0.0f, 0.0f, constants.g);

    // Heading reference. Always north for aborts (TODO: Maybe worth switching to actual roll in the future. for now not worth)
    Vector3 x_C(1.0f, 0.0f, 0.0f);

    // Attitude, angular velocity, angular acceleration, and specific force
    Vector3 e = P_acc + G_vec;
    Vector3 ed = P_jerk;
    Vector3 edd = P_snap;

    float n = e.norm();
    if (n < 1e-6f) {
        e = Vector3(0.0f, 0.0f, 1e-6f);
        n = 1e-6f;
    }
    Vector3 u = e / n;
    float nd = u.dot(ed);
    Vector3 ud = (ed - u * nd) / n;
    float ndd = ud.dot(ed) + u.dot(edd);
    Vector3 udd = (edd - 2.0f * ud * nd - u * ndd) / n;

    // Body frame TRIAD
    Vector3 x_ref = x_C;
    Vector3 c = u.cross(x_C);
    float nc = c.norm();
    if (nc < 1e-6f) {
        Vector3 x_C_alt(0.0f, 1.0f, 0.0f);
        c = u.cross(x_C_alt);
        nc = c.norm();
        x_ref = x_C_alt;
    }
    Vector3 b2 = c / nc;
    Vector3 b1 = b2.cross(u);

    // First and second time derivatives of b2 and b1
    Vector3 cd = ud.cross(x_ref);
    float ncd = b2.dot(cd);
    Vector3 b2d = (cd - b2 * ncd) / nc;
    Vector3 cdd = udd.cross(x_ref);
    float ncdd = b2d.dot(cd) + b2.dot(cdd);
    Vector3 b2dd = (cdd - 2.0f * b2d * ncd - b2 * ncdd) / nc;

    Vector3 b1d = b2d.cross(u) + b2.cross(ud);
    Vector3 b1dd = b2dd.cross(u) + 2.0f * b2d.cross(ud) + b2.cross(udd);

    // Body angular velocity w_B and angular acceleration wd_B
    float wx = -ud.dot(b2);
    float wy = ud.dot(b1);
    float wz = b1d.dot(b2);

    float wxd = -udd.dot(b2) - ud.dot(b2d);
    float wyd = udd.dot(b1) + ud.dot(b1d);
    float wzd = b1dd.dot(b2) + b1d.dot(b2d);

    Vector3 w_B(wx, wy, wz);
    Vector3 wd_B(wxd, wyd, wzd);

    // Quaternion from DCM C_IB = [b1, b2, b3]
    Matrix3_3 C_IB;
    C_IB << b1, b2, u;
    Vector4 qk = DCM_Quat_Conversion(C_IB);
    if (q_prev != nullptr && q_prev->dot(qk) < 0.0f) {
        qk = -qk;
    }

    // CG kinematics from P-space
    Vector3 r_cg = P_pos - rho * u;
    Vector3 v_cg = P_vel - rho * ud;
    AccelCg_out = P_acc - rho * udd;

    // Rotational Euler dynamics
    Vector3 J_w = constants.J.cwiseProduct(w_B);
    Vector3 tau_B = constants.J.cwiseProduct(wd_B) + w_B.cross(J_w);

    // Exact body force required to track accel
    Matrix3_3 C_BI = C_IB.transpose();
    Vector3 F_B = constants.m * (C_BI * (AccelCg_out + G_vec));

    X_out << qk, r_cg, v_cg, w_B, 0.0f, 0.0f;
    Wrench_out.F_body = F_B;
    Wrench_out.Tau_body = tau_B;
}
