// PolyQP/PolyQP_types.h first: it defines PolyQPConstants (ComputeJtot) and configures Eigen for PolyQP before matlab_funcs.h includes it.
#include "PolyQP/PolyQP_types.h"
#include "matlab_funcs.h"

Matrix3_3 zetaCross(Vector3 zeta) {
    Matrix3_3 zc;
    zc << 0, -zeta[2], zeta[1], zeta[2], 0, -zeta[0], -zeta[1], zeta[0], 0;
    return zc;
}

Matrix3_3 quatRot(Vector4 q) {
    Matrix3_3 CIB;
    CIB << 1 - 2 * (q[2] * q[2] + q[3] * q[3]), 2 * (q[1] * q[2] + q[0] * q[3]), 2 * (q[1] * q[3] - q[0] * q[2]), 2 * (q[1] * q[2] - q[0] * q[3]), 1 - 2 * (q[1] * q[1] + q[3] * q[3]),
    2 * (q[2] * q[3] + q[0] * q[1]), 2 * (q[1] * q[3] + q[0] * q[2]), 2 * (q[2] * q[3] - q[0] * q[1]), 1 - 2 * (q[1] * q[1] + q[2] * q[2]);
    return CIB;
}

Matrix4_4 HamiltonianProd(Vector4 q) {
    Matrix4_4 M;
    M << q(0), -q(1), -q(2), -q(3), q(1), q(0), -q(3), q(2), q(2), q(3), q(0), -q(1), q(3), -q(2), q(1), q(0);
    return M;
}

Vector16 StateAUG(Vector13 XKF, Vector3 G) {
    // Change Filter State to Controls State
    Vector16 X;
    X << XKF.segment<4>(0), XKF.segment<3>(4), XKF.segment<3>(7), G - XKF.segment<3>(10), XKF.segment<3>(10);
    return X;
}

// Convert a DCM to a quaternion robustly
Vector4 DCM_Quat_Conversion(Matrix3_3 R) {

    // Extract trace
    float tr = R(0, 0) + R(1, 1) + R(2, 2);

    // Pre-allocate output
    Vector4 q = Vector4::Zero();

    // CASE 1: Standard Case (Trace > 0)
    // This works for most orientations (small angles, < 90 deg tilt)
    if (tr > 0) {
        float S = sqrt(tr + 1.0) * 2; // S = 4 * qw
        q[0] = 0.25 * S;
        q[1] = (R(2, 1) - R(1, 2)) / S;
        q[2] = (R(0, 2) - R(2, 0)) / S;
        q[3] = (R(1, 0) - R(0, 1)) / S;
    }
    // CASE 2: Singularity Avoidance (Trace <= 0)
    // We must find the largest diagonal element to avoid dividing by zero.
    else {
        if ((R(0, 0) > R(1, 1)) && (R(0, 0) > R(2, 2))) {
            // Column 1 (X) is dominant
            float S = sqrt(1.0 + R(0, 0) - R(1, 1) - R(2, 2)) * 2;
            q[0] = (R(2, 1) - R(1, 2)) / S;
            q[1] = 0.25 * S;
            q[2] = (R(0, 1) + R(1, 0)) / S;
            q[3] = (R(0, 2) + R(2, 0)) / S;
        } else if (R(1, 1) > R(2, 2)) {
            // Column 2 (Y) is dominant
            float S = sqrt(1.0 + R(1, 1) - R(0, 0) - R(2, 2)) * 2;
            q[0] = (R(0, 2) - R(2, 0)) / S;
            q[1] = (R(0, 1) + R(1, 0)) / S;
            q[2] = 0.25 * S;
            q[3] = (R(1, 2) + R(2, 1)) / S;
        } else {
            // Column 3 (Z) is dominant
            float S = sqrt(1.0 + R(2, 2) - R(0, 0) - R(1, 1)) * 2;
            q[0] = (R(1, 0) - R(0, 1)) / S;
            q[1] = (R(0, 2) + R(2, 0)) / S;
            q[2] = (R(1, 2) + R(2, 1)) / S;
            q[3] = 0.25 * S;
        }
    }

    // Enforce normalization to correct numerical drift
    q.normalize();

    return q;
}

// Transpiled from TOAD/Helper/ComputeJtot.m @ 5ac6bd23 (double; moved from cpp/PolyQP/ComputeJtot.cpp)
void ComputeJtot(double m_lox, double m_ipa, const PolyQP::PolyQPConstants &constants6DoF, Matrix3_3d &J_tot, double &CGz) {

    const double m_dry = constants6DoF.m_dry;
    const double m = m_dry + m_lox + m_ipa;

    // Propellant fill height
    double OxFluidHeight;
    double FuFluidHeight;
    if (constants6DoF.Vehicle == 0) {
        OxFluidHeight = 0;
        FuFluidHeight = 0;
    } else {
        OxFluidHeight = (m_lox / constants6DoF.OxMass) * constants6DoF.OxHeight * 0.9;
        FuFluidHeight = (m_ipa / constants6DoF.FuMass) * constants6DoF.FuHeight * 0.9;
    }

    // Propellant inertias (about their own fluid-column centroid)
    double J_xx = 1.0 / 12.0 * m_lox * (3 * std::pow(constants6DoF.OxRadius, 2.0) + std::pow(OxFluidHeight, 2.0));
    double J_zz = 1.0 / 2.0 * m_lox * std::pow(constants6DoF.OxRadius, 2.0);
    Matrix3_3d J_lox = Vector3d(J_xx, J_xx, J_zz).asDiagonal();

    J_xx = 1.0 / 12.0 * m_ipa * (3 * std::pow(constants6DoF.FuRadius, 2.0) + std::pow(FuFluidHeight, 2.0));
    J_zz = 1.0 / 2.0 * m_ipa * std::pow(constants6DoF.FuRadius, 2.0);
    Matrix3_3d J_ipa = Vector3d(J_xx, J_xx, J_zz).asDiagonal();

    // Fluid fill location & instantaneous CG (engine-attachment frame)
    const double OxFluidLocation = constants6DoF.Ox_Z + OxFluidHeight / 2;
    const double FuFluidLocation = constants6DoF.Fu_Z + FuFluidHeight / 2;
    CGz = (m_dry * constants6DoF.rTB + m_lox * OxFluidLocation + m_ipa * FuFluidLocation) / m;

    // Parallel-axis shift of each component to the instantaneous CG
    const double d_dry = constants6DoF.rTB - CGz;
    const double d_lox = OxFluidLocation - CGz;
    const double d_ipa = FuFluidLocation - CGz;

    const Matrix3_3d diag_dry = Vector3d(std::pow(d_dry, 2.0), std::pow(d_dry, 2.0), 0.0).asDiagonal();
    const Matrix3_3d J_dry = constants6DoF.J + m_dry * diag_dry;
    const Matrix3_3d diag_lox = Vector3d(std::pow(d_lox, 2.0), std::pow(d_lox, 2.0), 0.0).asDiagonal();
    J_lox = J_lox + m_lox * diag_lox;
    const Matrix3_3d diag_ipa = Vector3d(std::pow(d_ipa, 2.0), std::pow(d_ipa, 2.0), 0.0).asDiagonal();
    J_ipa = J_ipa + m_ipa * diag_ipa;

    // Bring everything to the instantaneous CG
    J_tot = J_dry + J_lox + J_ipa;
}
