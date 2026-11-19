#include "matlab_funcs.h"
#include "simstruc.h"
#include <math.h>

const float mag_bias_tau = 60.0f;

Matrix15_15 FlightStateTransitionMat2(Vector3 accel, Vector3 gyro, Matrix3_3 R_b2i) {
    Matrix15_15 F = Matrix15_15::Zero();
    F.block<3, 3>(0, 0) = -zetaCross(gyro);
    F.block<3, 3>(0, 9) = -Matrix3_3::Identity();
    F.block<3, 3>(3, 6) = Matrix3_3::Identity();
    F.block<3, 3>(6, 0) = -R_b2i * zetaCross(accel);
    F.block<3, 3>(12, 12) = -(1.0f / mag_bias_tau) * Matrix3_3::Identity();
    return F;
}

Vector19 FlightEstimator2(Vector19 x_est, constantsASTRA_t constantsASTRA, Vector15 z, float dT, Matrix15_15 &P, bool new_mag_packet, bool new_gps_packet, Matrix3_3 gps_vel_covar, Matrix3_3 gps_pos_covar) {
    // M-EKF Implementation
    // Remove bias from IMU and magnetometer
    z.segment<3>(0) = z.segment<3>(0) - x_est.segment<3>(13);
    z.segment<3>(3) = z.segment<3>(3) - x_est.segment<3>(10);
    z.segment<3>(6) = z.segment<3>(6) - x_est.segment<3>(16);
    float S = fmaxf(z.segment<3>(6).norm(), 1e-4f);
    z.segment<3>(6) = z.segment<3>(6) / S;

    // Extract quaternion
    Vector15 dx = Vector15::Zero();
    Vector4 q = x_est.segment<4>(0);
    float gyro_norm = z.segment<3>(3).norm();
    if (gyro_norm < 1e-5f || dT == 0.0f) {
        x_est.segment<4>(0) = q;
    } else {
        float half_angle = gyro_norm * dT * 0.5f;
        Vector4 delta_q;
        delta_q << cosf(half_angle), (z.segment<3>(3) / gyro_norm) * sinf(half_angle);
        x_est.segment<4>(0) = HamiltonianProd(q) * delta_q;
    }
    q = x_est.segment<4>(0) / x_est.segment<4>(0).norm();

    // A-priori quaternion estimate and rotation matrix
    q.normalize();
    Matrix3_3 R_b2i = quatRot(q).transpose();

    // GPS velocity and position correction
    Vector3 rGPS = (Vector3() << 0.0f, 0.0f, 0.31f).finished();
    z.segment<3>(12) = z.segment<3>(12) - R_b2i * z.segment<3>(3).cross(rGPS);
    z.segment<3>(9) = z.segment<3>(9) - R_b2i * rGPS;

    // State Transition Matrix
    Matrix15_15 F = FlightStateTransitionMat2(z.segment<3>(0), z.segment<3>(3), R_b2i);

    // Propagate rest of state using IMU (and OU mag bias decay)
    x_est.segment<3>(7) = x_est.segment<3>(7) + (R_b2i * z.segment<3>(0) - (Vector3() << 0.0f, 0.0f, constantsASTRA.g).finished()) * dT;
    x_est.segment<3>(4) = x_est.segment<3>(4) + x_est.segment<3>(7) * dT;
    x_est.segment<3>(16) = x_est.segment<3>(16) * expf(-dT / mag_bias_tau);

    // Discrete STM
    Matrix15_15 Phi = matrixExpPade6<Matrix15_15>(F * dT);

    // Magnetometer normalization factor
    Matrix3_3 MagMatrix = (Matrix3_3::Identity() - z.segment<3>(6) * z.segment<3>(6).transpose()) / S;

    // Extract Matrices
    Matrix18_18 Q_full = Q_gen(dT);
    Matrix15_15 Q = Matrix15_15::Zero();
    Q.block<12, 12>(0, 0) = Q_full.block<12, 12>(0, 0);
    Q.block<3, 3>(12, 12) = Q_full.block<3, 3>(15, 15);

    // Unit-sphere angular mag covariance: (sigma_m / S)^2 + sigma_dist^2
    float sigma_m = 0.009f;
    float sigma_dist_rad = 0.06f;
    float ratio = sigma_m / S;
    float sigma_u2 = ratio * ratio + sigma_dist_rad * sigma_dist_rad;
    Matrix3_3 R_mag = sigma_u2 * Matrix3_3::Identity();

    // Process Noise Covariance and a-priori propagation step
    P = Phi * P * Phi.transpose() + Q;
    Matrix15_15 p_new = (P + P.transpose()) / 2.0f;
    P = p_new;

    // IMU Update
    if (new_mag_packet) {
        // Predicted measurements
        Vector3 z_hat = R_b2i.transpose() * constantsASTRA.mag.normalized();

        // Measurement matrix
        Matrix3_15 H = Matrix3_15::Zero();
        H.block<3, 3>(0, 0) = zetaCross(z_hat);
        H.block<3, 3>(0, 12) = MagMatrix;

        // A priori covariance and Kalman gain
        Matrix15_3 L = P * H.transpose() * (H * P * H.transpose() + R_mag).inverse();

        // Kalman Gain Weighting based on predicted magnetic field
        Matrix15_15 ILH = (Matrix15_15::Identity() - L * H);
        P = ILH * P * ILH.transpose() + L * R_mag * L.transpose();
        Matrix15_15 p_new_gain = (P + P.transpose()) / 2.0f;
        P = p_new_gain;
        Vector3 residual = (z.segment<3>(6) - z_hat);
        dx = dx + L * residual;
    }

    // GPS Update
    if (new_gps_packet) {
        // Measurement matrix
        Matrix6_15 H = Matrix6_15::Zero();
        H.block<3, 3>(0, 0) = -R_b2i * zetaCross(rGPS);
        H.block<3, 3>(0, 3) = Matrix3_3::Identity();
        H.block<3, 3>(3, 0) = -R_b2i * zetaCross(z.segment<3>(3).cross(rGPS));
        H.block<3, 3>(3, 6) = Matrix3_3::Identity();
        H.block<3, 3>(3, 9) = R_b2i * zetaCross(rGPS);

        // Measurement Covariance Matrix
        Matrix6_6 R_gps = Matrix6_6::Zero();
        R_gps.block<3, 3>(0, 0) = 0.01f * Matrix3_3::Identity();
        R_gps.block<3, 3>(3, 3) = 0.01f * Matrix3_3::Identity();

        // A priori covariance and Kalman gain
        Matrix15_6 L = P * H.transpose() * (H * P * H.transpose() + R_gps).inverse();

        // Predicted measurements
        Vector6 z_hat;
        z_hat << x_est.segment<3>(4), x_est.segment<3>(7);

        // Kalman Gain Weighting based on predicted GPS
        Matrix15_15 ILH = (Matrix15_15::Identity() - L * H);
        P = ILH * P * ILH.transpose() + L * R_gps * L.transpose();
        Matrix15_15 p_new_gps = (P + P.transpose()) / 2.0f;
        P = p_new_gps;
        Vector6 residual = (z.segment<6>(9) - z_hat) - H * dx;
        Vector15 inn = L * residual;
        dx = dx + inn;
    }

    // Update full-state estimates
    Vector4 dq = (Vector4() << 1.0f, dx.segment<3>(0) / 2.0f).finished();
    dq.normalize();

    Vector4 q_nom = HamiltonianProd(q) * dq;
    q_nom.normalize();
    x_est.segment<4>(0) = q_nom;
    x_est.segment<9>(4) = x_est.segment<9>(4) + dx.segment<9>(3);
    x_est.segment<3>(16) = x_est.segment<3>(16) + dx.segment<3>(12);
    return x_est;
}
