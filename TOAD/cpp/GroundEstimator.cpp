#include "matlab_funcs.h"
#define RTK 1

Matrix18_18 GroundStateTransitionMat(Vector3 accel, Vector3 gyro, Matrix3_3 R_b2i) {
  // Remove angular rates from error-state. Make safety copies of all relevant
  // files into Archive
  Matrix18_18 F = Matrix18_18::Zero();
  F.block<3, 3>(0, 0) = -zetaCross(gyro);
  F.block<3, 3>(0, 9) = -Matrix3_3::Identity();
  F.block<3, 3>(3, 6) = Matrix3_3::Identity();
  F.block<3, 3>(6, 0) = -R_b2i * zetaCross(accel);
  F.block<3, 3>(6, 12) = -R_b2i;
  return F;
}

Vector19 GroundEstimator(Vector19 x_est, constantsASTRA_t constantsASTRA, Vector15 z, float dT, Matrix18_18 &P, bool new_mag_packet, bool new_gps_packet, Matrix3_3 gps_vel_covar,
                         Matrix3_3 gps_pos_covar) {
  // M-EKF Implementation
  // Remove bias from IMU
  z.segment<3>(0) = z.segment<3>(0) - x_est.segment<3>(13);
  z.segment<3>(3) = z.segment<3>(3) - x_est.segment<3>(10);
  z.segment<3>(6) = z.segment<3>(6) - x_est.segment<3>(16);
  float mag_norm = z.segment<3>(6).norm();
  z.segment<3>(6).normalize();

  // Extract quaternion
  Vector18 dx = Vector18::Zero();
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

  // State Transition Matrix
  Matrix18_18 F = GroundStateTransitionMat(z.segment<3>(0), z.segment<3>(3), R_b2i);

  // Propagate rest of state using IMU
  x_est.segment<3>(7) = x_est.segment<3>(7) + (R_b2i * z.segment<3>(0) - (Vector3() << 0, 0, constantsASTRA.g).finished()) * dT;
  x_est.segment<3>(4) = x_est.segment<3>(4) + x_est.segment<3>(7) * dT;

  // Discrete STM
  Matrix18_18 Phi = matrixExpPade6<Matrix18_18>(F * dT);

  // Extract Matrices
  Matrix18_18 Q = constantsASTRA.Q;
  Matrix6_6 R = constantsASTRA.R;

  // Process Noise Covariance and a-priori propagation step
  // Q = 0.5 * Q;
  Matrix3_3 mag_matrix = (Matrix3_3::Identity() - z.segment<3>(6) * z.segment<3>(6).transpose()) / mag_norm;
  R.block<3, 3>(3, 3) = 4e-1 * mag_matrix + 5e-5 * z.segment<3>(6) * z.segment<3>(6).transpose();
  P = Phi * P * Phi.transpose() + Q;
  P = (P + P.transpose()) / 2;

  // IMU Update
  if (new_mag_packet) {
    // Measurement matrix
    Matrix6_18 H = Matrix6_18::Zero();
    H.block<3, 3>(0, 0) = zetaCross(R_b2i.transpose() * (Vector3() << 0, 0, constantsASTRA.g).finished());
    H.block<3, 3>(0, 12) = Matrix3_3::Identity();
    H.block<3, 3>(3, 0) = zetaCross(R_b2i.transpose() * constantsASTRA.mag);
    H.block<3, 3>(3, 15) = mag_matrix;

    // A priori covariance and Kalman gain
    Matrix18_6 L = P * H.transpose() * (H * P * H.transpose() + R).inverse();
    // L.block<1, 6>(14, 0) = Vector6::Zero().transpose();

    // Predicted measurements
    Vector6 z_hat = (Vector6() << R_b2i.transpose() * (Vector3() << 0, 0, constantsASTRA.g).finished(), R_b2i.transpose() * constantsASTRA.mag).finished();

    // Kalman Gain Weighting based on predicted acceleration
    Matrix18_18 ILH = (Matrix18_18::Identity() - L * H);
    P = ILH * P * ILH.transpose() + L * R * L.transpose();
    P = (P + P.transpose()) / 2;
    Vector6 residual = ((Vector6() << z.segment<3>(0), z.segment<3>(6)).finished() - z_hat);
    dx = dx + L * residual;
  }

  // GPS Update
  if (new_gps_packet) {
    // Measurement matrix
    Matrix6_18 H = Matrix6_18::Zero();
    H.block<3, 3>(0, 3) = Matrix3_3::Identity();
    H.block<3, 3>(3, 6) = Matrix3_3::Identity();

    Matrix6_6 R = Matrix6_6::Zero();
    R.block<3, 3>(0, 0) = 0.05 * Matrix3_3::Identity();
    R.block<3, 3>(3, 3) = 0.05 * Matrix3_3::Identity();

    // A priori covariance and Kalman gain
    Matrix18_6 L = P * H.transpose() * (H * P * H.transpose() + R).inverse();
    // L.block<1, 6>(14, 0) = Vector6::Zero().transpose();

    // Predicted measurements
    Vector6 z_hat = (Vector6() << x_est.segment<3>(4), x_est.segment<3>(7)).finished();

    // Kalman Gain Weighting based on predicted acceleration
    Matrix18_18 ILH = (Matrix18_18::Identity() - L * H);
    P = ILH * P * ILH.transpose() + L * R * L.transpose();
    P = (P + P.transpose()) / 2;
    Vector6 residual = (z.segment<6>(9) - z_hat);
    Vector18 inn = L * residual;
    dx = dx + inn;
  }

  // Update full-state estimates
  Vector4 dq = (Vector4() << 1, dx.segment<3>(0) / 2).finished();
  dq.normalize();

  Vector4 q_nom = HamiltonianProd(q) * dq;
  q_nom.normalize();
  x_est.segment<4>(0) = q_nom;
  x_est.segment<15>(4) = x_est.segment<15>(4) + dx.segment<15>(3);
  x_est.segment<6>(4) = Vector6::Zero();

  return x_est;
}