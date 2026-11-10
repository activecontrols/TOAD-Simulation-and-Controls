#include "matlab_funcs.h"

#define RTK 1

Matrix9_9 FlightStateTransitionMat(Vector3 accel, Vector3 gyro, Matrix3_3 R_b2i) {
    // Remove angular rates from error-state. Make safety copies of all relevant
    // files into Archive
    Matrix9_9 F = Matrix9_9::Zero();
    F.block<3, 3>(0, 0) = -zetaCross(gyro);
    F.block<3, 3>(3, 6) = Matrix3_3::Identity();
    F.block<3, 3>(6, 0) = -R_b2i * zetaCross(accel);
    return F;
}

Vector19 FlightEstimator(Vector19 x_est, constantsASTRA_t constantsASTRA, Vector15 z, float dT, Matrix9_9 &P, bool new_gps_packet, Matrix3_3 gps_vel_covar, Matrix3_3 gps_pos_covar) {
    // M-EKF Implementation
    // Remove bias from IMU
    z.segment<3>(0) = z.segment<3>(0) - x_est.segment<3>(13);
    z.segment<3>(3) = z.segment<3>(3) - x_est.segment<3>(10);
    z.segment<3>(6) = z.segment<3>(6) - x_est.segment<3>(16);

    // Extract quaternion
    Vector9 dx = Vector9::Zero();
    Vector4 q = x_est.segment<4>(0);
    Vector4 qdot = 0.5 * HamiltonianProd(q) * (Vector4() << 0, z.segment<3>(3)).finished();
    x_est.segment<4>(0) = q + qdot * dT;
    q = x_est.segment<4>(0);

    // A-priori quaternion estimate and rotation matrix
    q.normalize();
    Matrix3_3 R_b2i = quatRot(q).transpose();

    // GPS Gyroscopic Correction
    Vector3 rGPS = (Vector3() << 0, 0, 0.31).finished();
    z.segment<3>(12) = z.segment<3>(12) - R_b2i * z.segment<3>(3).cross(rGPS);
    z.segment<3>(9) = z.segment<3>(9) - R_b2i * rGPS;

    // State Transition Matrix
    Matrix9_9 F = FlightStateTransitionMat(z.segment<3>(0), z.segment<3>(3), R_b2i);

    // Propagate rest of state using IMU
    x_est.segment<3>(7) = x_est.segment<3>(7) + (R_b2i * z.segment<3>(0) - (Vector3() << 0, 0, constantsASTRA.g).finished()) * dT;
    x_est.segment<3>(4) = x_est.segment<3>(4) + x_est.segment<3>(7) * dT;

    // Discrete STM
    Matrix9_9 Phi = matrixExpPade6<Matrix9_9>(F * dT);

    // Extract Matrices
    Matrix9_9 Q = constantsASTRA.Q.block<9, 9>(0, 0);

    // Process Noise Covariance and a-priori propagation step
    Q = 0.5 * Q;
    P = Phi * P * Phi.transpose() + Q;

    // GPS Update
    if (new_gps_packet) {
        // Measurement matrix
        Matrix6_9 H = Matrix6_9::Zero();
        H.block<3, 3>(0, 3) = Matrix3_3::Identity();
        H.block<3, 3>(3, 6) = Matrix3_3::Identity();

        Matrix6_6 R = Matrix6_6::Zero();
        R.block<3, 3>(0, 0) = 1.5 * Matrix3_3::Identity(); // TODO - use covariances
        R.block<3, 3>(3, 3) = 0.1 * Matrix3_3::Identity();

        // A priori covariance and Kalman gain
        Matrix9_6 L = P * H.transpose() * (H * P * H.transpose() + R).inverse();

        // Predicted measurements
        Vector6 z_hat = (Vector6() << x_est.segment<3>(4), x_est.segment<3>(7)).finished();

        // Update Error State
        Matrix9_9 ILH = (Matrix9_9::Identity() - L * H);
        P = ILH * P * ILH.transpose() + L * R * L.transpose();
        Vector6 residual = (z.segment<6>(9) - z_hat);
        Vector9 inn = L * residual;
        dx = dx + inn;
    }

    // Update full-state estimates
    Vector4 dq = (Vector4() << 1, dx.segment<3>(0) / 2).finished();
    dq.normalize();

    Vector4 q_nom = HamiltonianProd(q) * dq;
    q_nom.normalize();
    x_est.segment<4>(0) = q_nom;
    x_est.segment<6>(4) = x_est.segment<6>(4) + dx.segment<6>(3);
    return x_est;
}