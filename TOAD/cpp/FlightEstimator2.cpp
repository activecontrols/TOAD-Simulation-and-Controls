#include "matlab_funcs.h"
#include "simstruc.h"
#include <cmath>

static constexpr float kMagBiasTau = 60.0f; // Ornstein-Uhlenbeck damping time constant (s)

static Matrix15_15 FlightStateTransitionMat(Vector3 accel, Vector3 gyro, Matrix3_3 R_b2i) {
    Matrix15_15 F = Matrix15_15::Zero();
    F.block<3, 3>(0, 0)   = -zetaCross(gyro);
    F.block<3, 3>(0, 9)   = -Matrix3_3::Identity();
    F.block<3, 3>(3, 6)   = Matrix3_3::Identity();
    F.block<3, 3>(6, 0)   = -R_b2i * zetaCross(accel);
    F.block<3, 3>(12, 12) = -(1.0f / kMagBiasTau) * Matrix3_3::Identity();
    return F;
}

Vector19 FlightEstimator2(Vector19 x_est, constantsASTRA_t constantsASTRA, Vector15 z, float dT, Matrix15_15 &P, bool new_mag_packet, bool new_gps_packet, Matrix3_3 gps_vel_covar,
                          Matrix3_3 gps_pos_covar) {
    (void)gps_vel_covar;
    (void)gps_pos_covar;

    // Remove bias from IMU & mag, then project mag onto unit sphere
    z.segment<3>(0) -= x_est.segment<3>(13);
    z.segment<3>(3) -= x_est.segment<3>(10);
    Vector3 mag_corr = z.segment<3>(6) - x_est.segment<3>(16);
    float S = fmaxf(mag_corr.norm(), 1e-4f);
    Vector3 u_mag = mag_corr / S;

    // Quaternion propagation
    Vector15 dx = Vector15::Zero();
    Vector4 q = x_est.segment<4>(0);
    Vector3 omega = z.segment<3>(3);
    float gyro_norm = omega.norm();
    if (gyro_norm >= 1e-5f && dT > 0.0f) {
        float half_angle = gyro_norm * dT * 0.5f;
        Vector4 delta_q;
        delta_q << cosf(half_angle), (omega / gyro_norm) * sinf(half_angle);
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
        Matrix15_15 F = FlightStateTransitionMat(z.segment<3>(0), omega, R_b2i);
        x_est.segment<3>(7)  += (R_b2i * z.segment<3>(0) - Vector3(0.0f, 0.0f, constantsASTRA.g)) * dT;
        x_est.segment<3>(4)  += x_est.segment<3>(7) * dT;
        x_est.segment<3>(16) *= expf(-dT / kMagBiasTau);

        Matrix15_15 Phi = matrixExpPade6<Matrix15_15>(F * dT);
        Matrix18_18 Q_full = Q_gen(dT);
        Matrix15_15 Q = Matrix15_15::Zero();
        Q.block<12, 12>(0, 0) = Q_full.block<12, 12>(0, 0);
        Q.block<3, 3>(12, 12) = Q_full.block<3, 3>(15, 15);

        P = Phi * P * Phi.transpose() + Q;
        P = 0.5f * (P + P.transpose());
    }

    // Magnetometer update (unit-sphere angular covariance + slow OU mag-bias slack)
    if (new_mag_packet) {
        Vector3 mag_ref = constantsASTRA.mag.normalized();
        Vector3 z_hat = R_b2i.transpose() * mag_ref;

        Matrix3_15 H = Matrix3_15::Zero();
        H.block<3, 3>(0, 0)  = zetaCross(z_hat);
        H.block<3, 3>(0, 12) = (Matrix3_3::Identity() - u_mag * u_mag.transpose()) / S;

        const float sigma_m = 0.009f;
        const float sigma_dist_rad = 0.06f;
        float ratio = sigma_m / S;
        float sigma_u2 = ratio * ratio + sigma_dist_rad * sigma_dist_rad;
        Matrix3_3 R_mag = sigma_u2 * Matrix3_3::Identity();

        Matrix15_3 L = P * H.transpose() * (H * P * H.transpose() + R_mag).inverse();
        Matrix15_15 ILH = Matrix15_15::Identity() - L * H;
        P = ILH * P * ILH.transpose() + L * R_mag * L.transpose();
        P = 0.5f * (P + P.transpose());

        dx += L * (u_mag - z_hat);
    }

    // GPS position & velocity update with attitude/gyro-bias lever-arm Jacobians
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
    x_est.segment<9>(4)  += dx.segment<9>(3);
    x_est.segment<3>(16) += dx.segment<3>(12);
    return x_est;
}
