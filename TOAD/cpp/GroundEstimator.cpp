#include "matlab_funcs.h"
#include "simstruc.h"
#include <math.h>

Matrix15_15 GroundStateTransitionMat(Vector3 accel, Vector3 gyro, Matrix3_3 R_b2i) {
  // Error state: [dtheta(0:2), dr(3:5), dv(6:8), dbg(9:11), dba(12:14)]
  Matrix15_15 F = Matrix15_15::Zero();
  F.block<3, 3>(0, 0) = -zetaCross(gyro);
  F.block<3, 3>(0, 9) = -Matrix3_3::Identity();
  F.block<3, 3>(3, 6) = Matrix3_3::Identity();
  F.block<3, 3>(6, 0) = -R_b2i * zetaCross(accel);
  F.block<3, 3>(6, 12) = -R_b2i;
  return F;
}

Vector19 GroundEstimator(Vector19 x_est, constantsASTRA_t constantsASTRA, Vector15 z, float dT, Matrix15_15 &P, bool new_mag_packet, bool new_gps_packet, Matrix3_3 gps_vel_covar,
                         Matrix3_3 gps_pos_covar) {
  // M-EKF Implementation
  // Remove bias from IMU (mag bias is held 0 on stationary pad)
  z.segment<3>(0) = z.segment<3>(0) - x_est.segment<3>(13);
  z.segment<3>(3) = z.segment<3>(3) - x_est.segment<3>(10);
  float mag_norm = fmaxf(z.segment<3>(6).norm(), 1e-4f);
  z.segment<3>(6) = z.segment<3>(6) / mag_norm;

  // Extract quaternion
  Vector15 dx = Vector15::Zero();
  Vector4 q = x_est.segment<4>(0);
  Vector3 eulerAngles = z.segment<3>(3);
  if ((eulerAngles.norm() < 0.00001f) || dT == 0.0f) {
    q = x_est.segment<4>(0);
  } else {
    x_est.segment<4>(0) = (HamiltonianProd(q) * (Vector4() << cosf(eulerAngles.norm() * dT / 2.0f), eulerAngles.segment<3>(0) / (eulerAngles.norm()) * sinf(eulerAngles.norm() * dT / 2.0f)).finished());
    q = x_est.segment<4>(0);
  }

  // A-priori quaternion estimate and rotation matrix
  q.normalize();
  Matrix3_3 R_b2i = quatRot(q).transpose();

  // GPS lever-arm correction
  Vector3 rGPS = (Vector3() << 0.0f, 0.0f, 0.31f).finished();
  z.segment<3>(9) = z.segment<3>(9) - R_b2i * rGPS;
  z.segment<3>(12) = z.segment<3>(12) - R_b2i * z.segment<3>(3).cross(rGPS);

  // State Transition Matrix
  Matrix15_15 F = GroundStateTransitionMat(z.segment<3>(0), z.segment<3>(3), R_b2i);

  // Propagate rest of state using IMU
  x_est.segment<3>(7) = x_est.segment<3>(7) + (R_b2i * z.segment<3>(0) - (Vector3() << 0.0f, 0.0f, constantsASTRA.g).finished()) * dT;
  x_est.segment<3>(4) = x_est.segment<3>(4) + x_est.segment<3>(7) * dT;

  // Discrete STM
  Matrix15_15 Phi = matrixExpPade6<Matrix15_15>(F * dT);

  // Extract Matrices
  Matrix15_15 Q = Q_gen(dT).block<15, 15>(0, 0);

  // Process Noise Covariance and a-priori propagation step
  P = Phi * P * Phi.transpose() + Q;

  Matrix15_15 p_new = (P + P.transpose()) / 2.0f;
  P = p_new;

  // IMU Update (accel + yaw-only mag decoupled from roll/pitch)
  if (new_mag_packet) {
    Vector3 g_B = R_b2i.transpose() * (Vector3() << 0.0f, 0.0f, constantsASTRA.g).finished();
    Vector3 m_I_meas = R_b2i * z.segment<3>(6);
    float m_horiz = fmaxf(constantsASTRA.mag.head<2>().norm(), 1e-4f);
    float y_psi = atan2f(m_I_meas(0) * constantsASTRA.mag(1) - m_I_meas(1) * constantsASTRA.mag(0),
                         m_I_meas(0) * constantsASTRA.mag(0) + m_I_meas(1) * constantsASTRA.mag(1));

    // Measurement matrix
    Matrix4_15 H = Matrix4_15::Zero();
    H.block<3, 3>(0, 0) = zetaCross(g_B);
    H.block<3, 3>(0, 12) = Matrix3_3::Identity();
    H.block<1, 3>(3, 0) = R_b2i.row(2);

    float sigma_m = 0.009f;
    float sigma_dist_rad = 0.03f;
    float ratio = sigma_m / mag_norm;
    float sigma_u2 = ratio * ratio + sigma_dist_rad * sigma_dist_rad;
    Matrix4_4 R = Matrix4_4::Zero();
    R.block<3, 3>(0, 0) = 0.002f * Matrix3_3::Identity();
    R(3, 3) = sigma_u2 / (m_horiz * m_horiz);

    // A priori covariance and Kalman gain
    Matrix15_4 L = P * H.transpose() * (H * P * H.transpose() + R).inverse();

    // Predicted measurements
    Vector4 residual = (Vector4() << (z.segment<3>(0) - g_B), y_psi).finished();

    // Kalman Gain Weighting based on predicted acceleration
    Matrix15_15 ILH = (Matrix15_15::Identity() - L * H);
    P = ILH * P * ILH.transpose() + L * R * L.transpose();
    Matrix15_15 p_new_imu = (P + P.transpose()) / 2.0f;
    P = p_new_imu;
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

    Matrix6_6 R = Matrix6_6::Zero();
    R.block<3, 3>(0, 0) = 0.005f * Matrix3_3::Identity();
    R.block<3, 3>(3, 3) = 0.002f * Matrix3_3::Identity();

    // A priori covariance and Kalman gain
    Matrix15_6 L = P * H.transpose() * (H * P * H.transpose() + R).inverse();

    // Predicted measurements
    Vector6 z_hat = (Vector6() << x_est.segment<3>(4), x_est.segment<3>(7)).finished();

    // Kalman Gain Weighting based on predicted acceleration
    Matrix15_15 ILH = (Matrix15_15::Identity() - L * H);
    P = ILH * P * ILH.transpose() + L * R * L.transpose();
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
  x_est.segment<12>(4) = x_est.segment<12>(4) + dx.segment<12>(3);

  return x_est;
}
