#include "matlab_funcs.h"
#include "simstruc.h"
#include <cmath>

Matrix15_15 GroundStateTransitionMat(Vector3 accel, Vector3 gyro, Matrix3_3 R_b2i) {
    Matrix15_15 F = Matrix15_15::Zero();
    F.block<3, 3>(0, 0)  = -zetaCross(gyro);
    F.block<3, 3>(0, 9)  = -Matrix3_3::Identity();
    F.block<3, 3>(3, 6)  = Matrix3_3::Identity();
    F.block<3, 3>(6, 0)  = -R_b2i * zetaCross(accel);
    F.block<3, 3>(6, 12) = -R_b2i;
    return F;
}

Vector19 GroundEstimator(Vector19 x_est, constantsASTRA_t constantsASTRA, Vector15 z, float dT, Matrix15_15 &P, bool new_mag_packet, bool new_gps_packet, Matrix3_3 gps_vel_covar,
                         Matrix3_3 gps_pos_covar) {
    (void)gps_vel_covar;
    (void)gps_pos_covar;

    // Remove estimated IMU biases (mag bias is not estimated on stationary pad)
    z.segment<3>(0) -= x_est.segment<3>(13);
    z.segment<3>(3) -= x_est.segment<3>(10);
    float mag_norm = fmaxf(z.segment<3>(6).norm(), 1e-4f);
    Vector3 u_mag = z.segment<3>(6) / mag_norm;

    // Quaternion propagation
    Vector15 dx = Vector15::Zero();
    Vector4 q = x_est.segment<4>(0);
    Vector3 omega = z.segment<3>(3);
    float omega_norm = omega.norm();
    if (omega_norm >= 1e-5f && dT > 0.0f) {
        float half_angle = 0.5f * omega_norm * dT;
        Vector4 delta_q;
        delta_q << cosf(half_angle), (omega / omega_norm) * sinf(half_angle);
        q = HamiltonianProd(q) * delta_q;
    }
    q.normalize();
    Matrix3_3 R_b2i = quatRot(q).transpose();

    // GPS lever-arm compensation (position and velocity)
    Vector3 rGPS(0.0f, 0.0f, 0.31f);
    z.segment<3>(9)  -= R_b2i * rGPS;
    z.segment<3>(12) -= R_b2i * omega.cross(rGPS);

    // State and covariance propagation
    if (dT > 0.0f) {
        Matrix15_15 F = GroundStateTransitionMat(z.segment<3>(0), omega, R_b2i);
        x_est.segment<3>(7) += (R_b2i * z.segment<3>(0) - Vector3(0.0f, 0.0f, constantsASTRA.g)) * dT;
        x_est.segment<3>(4) += x_est.segment<3>(7) * dT;

        Matrix15_15 Phi = matrixExpPade6<Matrix15_15>(F * dT);
        Matrix15_15 Q = Q_gen(dT).block<15, 15>(0, 0);
        P = Phi * P * Phi.transpose() + Q;
        P = 0.5f * (P + P.transpose());
    }

    // Synchronized IMU Accel + Yaw-Only Magnetometer Update
    if (new_mag_packet) {
        Vector3 g_B = R_b2i.transpose() * Vector3(0.0f, 0.0f, constantsASTRA.g);

        Vector2 m_ref_xy(constantsASTRA.mag(0), constantsASTRA.mag(1));
        float m_horiz = fmaxf(m_ref_xy.norm(), 1e-4f);
        Vector2 b_ref = m_ref_xy / m_horiz;

        Vector3 m_I_meas = R_b2i * u_mag;
        Vector2 m_meas_xy(m_I_meas(0), m_I_meas(1));
        float m_meas_horiz = m_meas_xy.norm();

        // First-principles unit-sphere mag variance: (sigma_m / S)^2 + sigma_dist^2
        const float sigma_m = 0.009f;
        const float sigma_dist_rad = 0.03f;
        float ratio = sigma_m / mag_norm;
        float sigma_u2 = ratio * ratio + sigma_dist_rad * sigma_dist_rad;
        const float r_accel = 0.002f;

        if (m_meas_horiz > 1e-3f) {
            Vector2 b_meas = m_meas_xy / m_meas_horiz;
            float y_psi = atan2f(b_meas(0) * b_ref(1) - b_meas(1) * b_ref(0),
                                 b_meas(0) * b_ref(0) + b_meas(1) * b_ref(1));

            Matrix4_15 H = Matrix4_15::Zero();
            H.block<3, 3>(0, 0)  = zetaCross(g_B);
            H.block<3, 3>(0, 12) = Matrix3_3::Identity();
            H.block<1, 3>(3, 0)  = R_b2i.row(2); // Yaw-only: strictly decoupled from roll/pitch

            Matrix4_4 R_imu = Matrix4_4::Zero();
            R_imu.block<3, 3>(0, 0) = r_accel * Matrix3_3::Identity();
            R_imu(3, 3) = sigma_u2 / (m_horiz * m_horiz);

            Matrix15_4 L = P * H.transpose() * (H * P * H.transpose() + R_imu).inverse();
            Matrix15_15 ILH = Matrix15_15::Identity() - L * H;
            P = ILH * P * ILH.transpose() + L * R_imu * L.transpose();
            P = 0.5f * (P + P.transpose());

            Vector4 residual;
            residual << (z.segment<3>(0) - g_B), y_psi;
            dx += L * residual;
        } else {
            Matrix3_15 H = Matrix3_15::Zero();
            H.block<3, 3>(0, 0)  = zetaCross(g_B);
            H.block<3, 3>(0, 12) = Matrix3_3::Identity();

            Matrix3_3 R_acc = r_accel * Matrix3_3::Identity();
            Matrix15_3 L = P * H.transpose() * (H * P * H.transpose() + R_acc).inverse();
            Matrix15_15 ILH = Matrix15_15::Identity() - L * H;
            P = ILH * P * ILH.transpose() + L * R_acc * L.transpose();
            P = 0.5f * (P + P.transpose());

            dx += L * (z.segment<3>(0) - g_B);
        }
    }

    // GPS Position & Velocity Update
    if (new_gps_packet) {
        Matrix6_15 H = Matrix6_15::Zero();
        H.block<3, 3>(0, 0) = -R_b2i * zetaCross(rGPS);
        H.block<3, 3>(0, 3) = Matrix3_3::Identity();
        H.block<3, 3>(3, 0) = -R_b2i * zetaCross(omega.cross(rGPS));
        H.block<3, 3>(3, 6) = Matrix3_3::Identity();
        H.block<3, 3>(3, 9) = R_b2i * zetaCross(rGPS);

        Matrix6_6 R_gps = Matrix6_6::Zero();
        R_gps.block<3, 3>(0, 0) = 0.005f * Matrix3_3::Identity();
        R_gps.block<3, 3>(3, 3) = 0.002f * Matrix3_3::Identity();

        Matrix15_6 L = P * H.transpose() * (H * P * H.transpose() + R_gps).inverse();
        Matrix15_15 ILH = Matrix15_15::Identity() - L * H;
        P = ILH * P * ILH.transpose() + L * R_gps * L.transpose();
        P = 0.5f * (P + P.transpose());

        Vector6 z_hat;
        z_hat << x_est.segment<3>(4), x_est.segment<3>(7);
        Vector6 residual = (z.segment<6>(9) - z_hat) - H * dx;
        dx += L * residual;
    }

    // Inject error state into nominal state
    Vector4 dq;
    dq << 1.0f, 0.5f * dx.segment<3>(0);
    dq.normalize();

    Vector4 q_nom = HamiltonianProd(q) * dq;
    q_nom.normalize();
    x_est.segment<4>(0)   = q_nom;
    x_est.segment<12>(4) += dx.segment<12>(3);

    return x_est;
}
