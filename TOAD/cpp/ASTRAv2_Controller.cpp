#include "controller_and_estimator.h"
#include "matlab_funcs.h"
#include <cmath>

// LESO Augmented Time Varying LQR formulation
// New controlzzz

static Vector3 I_trans = Vector3::Zero();

void ASTRAv2_Controller_reset() {
  I_trans.setZero();
}

Vector4 ASTRAv2_Controller(bool GND, const Vector15 &X_est, const Vector15 &X_trg, const Vector4 &U_ff, const constantsASTRA_t &constantsASTRA, const Vector6 &Dist_LESO, float dT,
                           Controller_Internals *cs) {

  // Reset persistent states when on ground or first initialization
  if (dT <= 0.0f || GND == 1) {
    I_trans.setZero();
    return Vector4::Zero();
  }

  // Unpack gains and disturbance estimates
  Vector3 a_dist = Dist_LESO.segment<3>(0);
  Vector3 ang_accel_dist = Dist_LESO.segment<3>(3);

  float Mass = constantsASTRA.m + X_est(13) + X_est(14);
  Vector3 g_vec = (Vector3() << 0.0f, 0.0f, -constantsASTRA.g).finished();

  // Translational Trim (Outer Loop w / Integrator)
  Vector3 pos_err = X_est.segment<3>(4) - X_trg.segment<3>(4);
  Vector3 vel_err = X_est.segment<3>(7) - X_trg.segment<3>(7);

  // Leaked integration for the translational integral trim
  I_trans = I_trans * 0.999f + pos_err * dT;
  Vector3 MaxInt = (Vector3() << 0.35f, 0.35f, 0.15f).finished();
  I_trans = I_trans.cwiseMin(MaxInt).cwiseMax(-MaxInt);

  Vector6 trans_err;
  trans_err << pos_err, vel_err;
  Vector3 Delta_A = -constantsASTRA.K_trans * trans_err - constantsASTRA.K_I * I_trans;

  // Acceleration saturation limit
  Vector3 MaxAccelCorr = (Vector3() << 1.5f, 1.5f, 4.0f).finished();
  Delta_A = Delta_A.cwiseMin(MaxAccelCorr).cwiseMax(-MaxAccelCorr);

  // Nominal Trajectory Acceleration
  Vector4 Q_ref = X_trg.segment<4>(0);
  Vector4 Q_est = X_est.segment<4>(0);
  Matrix3_3 C_B2I_ref = quatRot(Q_ref).transpose();
  Matrix3_3 C_B2I_est = quatRot(Q_est).transpose();

  float theta_ff = U_ff[0];
  float phi_ff = U_ff[1];
  float thrust_ff = U_ff[2];

  Vector3 T_B_ff = thrust_ff * (Vector3() << cosf(theta_ff) * sinf(phi_ff), -sinf(theta_ff), cosf(theta_ff) * cosf(phi_ff)).finished();
  Vector3 a_ff = (C_B2I_ref * T_B_ff) / Mass + g_vec;

  // Commanded inertial acceleration with LESO disturbance compensation
  Vector3 a_cmd = a_ff - a_dist + Delta_A;
  Vector3 f_req = a_cmd - g_vec;

  // Dynamic Pullout Protection & Tilt Clamping
  // Protects and pulls the vehicle out of a trajectory if near unrecoverable
  float tb_norm = T_B_ff.norm();
  Vector3 Z_B_trg = C_B2I_ref * (T_B_ff / fmaxf(tb_norm, 1e-4f));
  float z_curr = X_est[6];
  float vz_curr = X_est[9];

  // Altitude & Sink - Rate Adaptive Cone Clamping
  // TODO (Pablo): Check if this is even worth having, or if u can simplify the logic on this.
  // If current velocity is negative indicating descent, we calculate the
  // stopping distance given our current altitude target and limit the
  // tilt accordingly.
  float MaxConeAng = 20.0f * M_PI / 180.0f;
  if (vz_curr < -1.0f) {
    float thrustDist = fminf(a_dist(2), 0.0f);
    float availAccel = fmaxf((constantsASTRA.MaxThrust / Mass + thrustDist) - constantsASTRA.g, 0.4f);
    float stopDist = 0.5f * (vz_curr * vz_curr) / availAccel;

    if (z_curr < stopDist * 1.35f + 2.5f) {
      // Emergency Braking Pullout: Lock thrust strictly vertical to maximize deceleration
      MaxConeAng = 10.0f * M_PI / 180.0f;
      f_req(2) = fmaxf(f_req(2), constantsASTRA.MaxThrust / Mass);
    } else if (z_curr < stopDist * 1.50f + 5.0f) {
      // Anticipatory Braking, squeeze cone angle smoothly
      float coneRatio = (z_curr - stopDist) / fmaxf(stopDist, 1.0f);
      float clampedRatio = fminf(fmaxf(coneRatio, 0.3f), 1.0f);
      MaxConeAng = (15.0f * M_PI / 180.0f) * clampedRatio;
      f_req(2) = fmaxf(f_req(2), constantsASTRA.g + 1.2f);
    } else {
      MaxConeAng = 20.0f * M_PI / 180.0f;
    }
  } else {
    MaxConeAng = 20.0f * M_PI / 180.0f;
  }

  // Project the required force direction into the allowable cone relative
  // to nominal thrust axis
  float norm_f = f_req.norm();
  float dotVal = f_req.dot(Z_B_trg) / fmaxf(norm_f * Z_B_trg.norm(), 1e-6f);
  if (dotVal < cosf(MaxConeAng)) {
    Vector3 f_perp = f_req - f_req.dot(Z_B_trg) * Z_B_trg;
    float norm_perp = f_perp.norm();
    if (norm_perp > 1e-6f) {
      f_perp /= norm_perp;
      f_req = (cosf(MaxConeAng) * Z_B_trg + sinf(MaxConeAng) * f_perp) * norm_f;
    } else {
      f_req = Z_B_trg * norm_f;
    }
    norm_f = f_req.norm();
  }

  Vector3 Z_t;
  if (norm_f > 1e-6f) {
    Z_t = f_req / norm_f;
  } else {
    Z_t = (Vector3() << 0.0f, 0.0f, 1.0f).finished();
  }

  // Actuator Disturbance Trim & Intermediate Frame C_T2B
  Vector3 torque = ang_accel_dist.cwiseProduct(constantsASTRA.J);
  Vector4 U_cmd = Vector4::Zero();

  // Correct input clamping to actual limits
  U_cmd[2] = fminf(fmaxf(norm_f * Mass, 0.2f * constantsASTRA.MaxThrust), constantsASTRA.MaxThrust);
  float u_cmd3 = fmaxf(U_cmd[2], 1e-3f);

  // Gimbal angles required to cancel CoM offset torque:
  Vector4 U_dist = Vector4::Zero();
  U_dist[0] = -1.0f * torque[0] / (constantsASTRA.rTB * u_cmd3);

  float cos_phi_approx = 1.0f - 0.5f * theta_ff * theta_ff;
  U_dist[1] = -1.0f * torque[1] / (constantsASTRA.rTB * cos_phi_approx * u_cmd3);
  U_dist[2] = torque[2];

  // Construct Thrust-to-Body Rotation Matrix (C_T2B)
  float theta_eff = theta_ff - U_dist(0);
  float phi_eff = phi_ff - U_dist(1);

  Matrix3_3 R_x, R_y;
  R_x << 1.0f, 0.0f, 0.0f, 0.0f, cosf(theta_eff), -sinf(theta_eff), 0.0f, sinf(theta_eff), cosf(theta_eff);

  R_y << cosf(phi_eff), 0.0f, sinf(phi_eff), 0.0f, 1.0f, 0.0f, -sinf(phi_eff), 0.0f, cosf(phi_eff);

  Matrix3_3 C_T2B = R_y * R_x;

  // Triad Gen
  Vector3 X_t_ref = C_B2I_ref * C_T2B.col(0);
  X_t_ref = X_t_ref - X_t_ref.dot(Z_t) * Z_t;

  Vector3 X_t, Y_t;
  if (X_t_ref.norm() > 1e-6f) {
    X_t_ref.normalize();
    Y_t = Z_t.cross(X_t_ref);
    Y_t.normalize();
    X_t = Y_t.cross(Z_t);
  } else {
    Y_t = (Vector3() << 0.0f, 1.0f, 0.0f).finished();
    X_t = Y_t.cross(Z_t);
    X_t.normalize();
    Y_t = Z_t.cross(X_t);
  }

  Matrix3_3 C_T2I_cmd;
  C_T2I_cmd << X_t, Y_t, Z_t;
  Matrix3_3 C_B2I_cmd = C_T2I_cmd * C_T2B.transpose();
  Vector4 Q_cmd = DCM_Quat_Conversion(C_B2I_cmd);
  Q_cmd.normalize();
  if (Q_cmd(0) < 0.0f) {
    Q_cmd = -Q_cmd;
  }

  // Rotational Inner-Loop Control
  Vector4 Q_cmd_Conj = (Vector4() << Q_cmd(0), -Q_cmd.segment<3>(1)).finished();
  Vector4 AttError = HamiltonianProd(Q_cmd_Conj) * X_est.segment<4>(0);
  if (AttError(0) < 0.0f) {
    AttError = -AttError;
  }

  Vector3 omegaTRG = C_B2I_est.transpose() * C_B2I_ref * X_trg.segment<3>(10);

  Vector3 att_vec = AttError.segment<3>(1);
  float n = att_vec.norm();
  Vector3 Axis = (n > 1e-6f) ? (att_vec / att_vec.norm()) : Vector3::Zero();
  Vector3 ThetaErr = 2.0f * atan2f(att_vec.norm(), AttError(0)) * Axis;

  Vector6 X_err_rot;
  X_err_rot << ThetaErr, (X_est.segment<3>(10) - omegaTRG);
  Vector3 Delta_u = -constantsASTRA.K_rot * X_err_rot;

  //  Actuator Trim Integration & Disturbance Allocation
  float thrustMax = constantsASTRA.MaxThrust;
  float gimbalMax = 10.0f * M_PI / 180.0f;

  // Feedback Trim Bounds
  Vector3 MaxTrim = (Vector3() << 10.0f * M_PI / 180.0f, 10.0f * M_PI / 180.0f, 3.0f).finished();
  Vector3 Trim_FB = Delta_u.cwiseMin(MaxTrim).cwiseMax(-MaxTrim);

  // Active Gimbal Disturbance Trim Feedforward
  Vector3 U_dist_chan = (Vector3() << U_dist(0), U_dist(1), U_dist(2)).finished();
  Vector3 U_dist_req = (-U_dist_chan).cwiseMin(MaxTrim).cwiseMax(-MaxTrim);

  // Combine feedback and disturbance trim
  Vector3 Trim = Trim_FB + U_dist_req;

  // Allocate to actuator channels
  Vector3 U_ff_chan = (Vector3() << U_ff(0), U_ff(1), U_ff(3)).finished();
  Vector3 Final_Channels = Trim + U_ff_chan;

  U_cmd(0) = fminf(fmaxf(Final_Channels(0), -gimbalMax), gimbalMax); // Pitch (rad)
  U_cmd(1) = fminf(fmaxf(Final_Channels(1), -gimbalMax), gimbalMax); // Yaw (rad)
  U_cmd(2) = fminf(fmaxf(U_cmd(2), 0.2f * thrustMax), thrustMax);    // Thrust (N)
  U_cmd(3) = fminf(fmaxf(Final_Channels(2), -4.0f), 4.0f);           // Roll (rad/s^2)

  return U_cmd;
}