#include "matlab_funcs.h"

Matrix12_12 FlightStateTransitionMat2(Vector3 accel, Vector3 gyro, Matrix3_3 R_b2i) {
    Matrix12_12 F = Matrix12_12::Zero();
    F.block<3, 3>(0, 0) = -zetaCross(gyro);
    F.block<3, 3>(0, 9) = -Matrix3_3::Identity();
    F.block<3, 3>(3, 6) = Matrix3_3::Identity();
    F.block<3, 3>(6, 0) = -R_b2i * zetaCross(accel);
    return F;
}

Vector19 FlightEstimator2(Vector19 x_est, constantsASTRA_t constantsASTRA, Vector15 z, float dT, Matrix12_12 &P, bool new_mag_packet, bool new_gps_packet, Matrix3_3 gps_vel_covar, Matrix3_3 gps_pos_covar) {
    // M-EKF Implementation
    // Remove bias from IMU
    z.segment<3>(0) = z.segment<3>(0) - x_est.segment<3>(13);
    z.segment<3>(3) = z.segment<3>(3) - x_est.segment<3>(10);
    z.segment<3>(6) = z.segment<3>(6) - x_est.segment<3>(16);
    float S = z.segment<3>(6).norm();
    z.segment<3>(6) = z.segment<3>(6) / S;

    // Extract quaternion
    Vector12 dx = Vector12::Zero();
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

    // GPS velocity correction
    Vector3 rGPS = (Vector3() << 0.0f, 0.0f, 0.31f).finished();
    z.segment<3>(12) = z.segment<3>(12) - R_b2i * z.segment<3>(3).cross(rGPS);

    // State Transition Matrix
    Matrix12_12 F = FlightStateTransitionMat2(z.segment<3>(0), z.segment<3>(3), R_b2i);

    // Propagate rest of state using IMU
    x_est.segment<3>(7) = x_est.segment<3>(7) + (R_b2i * z.segment<3>(0) - (Vector3() << 0.0f, 0.0f, constantsASTRA.g).finished()) * dT;
    x_est.segment<3>(4) = x_est.segment<3>(4) + x_est.segment<3>(7) * dT;

    // Discrete STM
    Matrix12_12 Phi = matrixExpPade6<Matrix12_12>(F * dT);

    // Magnetometer normalization factor
    Matrix3_3 MagMatrix = (Matrix3_3::Identity() - z.segment<3>(6) * z.segment<3>(6).transpose()) / S;

    // Extract Matrices
    Matrix12_12 Q = constantsASTRA.Q.block<12, 12>(0, 0);
    Matrix3_3 R_mag = constantsASTRA.R.block<3, 3>(3, 3);

    // Matrix updates for the mags
    R_mag = 1e-2f * MagMatrix + 1e-8f * (z.segment<3>(6) * z.segment<3>(6).transpose());

    // Process Noise Covariance and a-priori propagation step
    P = Phi * P * Phi.transpose() + Q;
    P = (P + P.transpose()) / 2.0f;

    // IMU Update
    if (new_mag_packet) {
        // Measurement matrix
        Matrix3_12 H = Matrix3_12::Zero();
        H.block<3, 3>(0, 0) = zetaCross(R_b2i.transpose() * constantsASTRA.mag);

        // Predicted measurements
        Vector3 z_hat = R_b2i.transpose() * constantsASTRA.mag;

        // A priori covariance and Kalman gain
        Matrix12_3 L = P * H.transpose() * (H * P * H.transpose() + R_mag).inverse();

        // Kalman Gain Weighting based on predicted acceleration
        Matrix12_12 ILH = (Matrix12_12::Identity() - L * H);
        P = ILH * P * ILH.transpose() + L * R_mag * L.transpose();
        P = (P + P.transpose()) / 2.0f;
        Vector3 residual = (z.segment<3>(6) - z_hat);
        dx = dx + L * residual;
    }

    // GPS Update
    if (new_gps_packet) {
        // Measurement matrix
        Matrix6_12 H = Matrix6_12::Zero();
        H.block<3, 3>(0, 3) = Matrix3_3::Identity();
        H.block<3, 3>(3, 6) = Matrix3_3::Identity();

        // Measurement Covariance Matrix
        Matrix3_3 GyroCovar = 1e-2f * Matrix3_3::Identity();
        Matrix6_6 R_gps = Matrix6_6::Zero();
        R_gps.block<3, 3>(0, 0) = 0.1f * Matrix3_3::Identity();
        Matrix3_3 cross_term = R_b2i * zetaCross(rGPS);
        R_gps.block<3, 3>(3, 3) = 0.3f * Matrix3_3::Identity() + cross_term * GyroCovar * cross_term.transpose();

        // A priori covariance and Kalman gain
        Matrix12_6 L = P * H.transpose() * (H * P * H.transpose() + R_gps).inverse();

        // Predicted measurements
        Vector6 z_hat;
        z_hat << x_est.segment<3>(4), x_est.segment<3>(7);

        // Kalman Gain Weighting based on predicted acceleration
        Matrix12_12 ILH = (Matrix12_12::Identity() - L * H);
        P = ILH * P * ILH.transpose() + L * R_gps * L.transpose();
        P = (P + P.transpose()) / 2.0f;
        Vector6 residual = (z.segment<6>(9) - z_hat) - H * dx;
        Vector12 inn = L * residual;
        dx = dx + inn;
    }

    // Update full-state estimates
    Vector4 dq = (Vector4() << 1.0f, dx.segment<3>(0) / 2.0f).finished();
    dq.normalize();

    Vector4 q_nom = HamiltonianProd(q) * dq;
    q_nom.normalize();
    x_est.segment<4>(0) = q_nom;
    x_est.segment<9>(4) = x_est.segment<9>(4) + dx.segment<9>(3);
    return x_est;
}
