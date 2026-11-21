#include "controller_and_estimator.h"
#include "matlab_funcs.h"
#include "simstruc.h"

namespace ControllerAndEstimator {
constantsASTRA_t constantsASTRA;
Matrix15_15 P;
Matrix15_15 Flight_P;
Vector19 x_est;
Vector2 lastEMA;
Matrix9_4 dnf_X;
Matrix9_4 dnf_Y;
Vector9 xhat_pos;
Vector6 xhat_att;
Vector4 last_U_cmd;
float last_thrust;
bool last_GND;
const float tau = 0.002; // seconds - time constant for ema low pass filter applied to gimbal angles
Vector15 last_X_ctrl = Vector15::Zero(); // abort harness: read-only copy of X_ctrl (get_controller_x_est)
Vector3 last_a_dist = Vector3::Zero();   // abort harness: read-only copy of a_dist (get_position_disturbance)

// Hand off ground attitude/pos/vel/gyro-bias covariance and reset OU mag-bias slack for flight.
static void init_flight_covariance() {
  Flight_P = Matrix15_15::Zero();
  Flight_P.block<12, 12>(0, 0) = P.block<12, 12>(0, 0);
  Flight_P.block<3, 3>(12, 12) = Matrix3_3::Identity() * (0.015f * 0.015f);
  x_est.segment<3>(16).setZero();
}

Matrix3_3 get_pos_cov_mtx(GPS_State gps_state) {
  Matrix3_3 pos_cov_mtx;
  pos_cov_mtx(0, 0) = gps_state.posCovNN;
  pos_cov_mtx(0, 1) = -gps_state.posCovNE;
  pos_cov_mtx(0, 2) = -gps_state.posCovND;
  pos_cov_mtx(1, 0) = -gps_state.posCovNE;
  pos_cov_mtx(1, 1) = gps_state.posCovEE;
  pos_cov_mtx(1, 2) = gps_state.posCovED;
  pos_cov_mtx(2, 0) = -gps_state.posCovND;
  pos_cov_mtx(2, 1) = gps_state.posCovED;
  pos_cov_mtx(2, 2) = gps_state.posCovDD;
  return pos_cov_mtx;
}

Matrix3_3 get_vel_cov_mtx(GPS_State gps_state) {
  Matrix3_3 vel_cov_mtx;
  vel_cov_mtx(0, 0) = gps_state.velCovNN;
  vel_cov_mtx(0, 1) = -gps_state.velCovNE;
  vel_cov_mtx(0, 2) = -gps_state.velCovND;
  vel_cov_mtx(1, 0) = -gps_state.velCovNE;
  vel_cov_mtx(1, 1) = gps_state.velCovEE;
  vel_cov_mtx(1, 2) = gps_state.velCovED;
  vel_cov_mtx(2, 0) = -gps_state.velCovND;
  vel_cov_mtx(2, 1) = gps_state.velCovED;
  vel_cov_mtx(2, 2) = gps_state.velCovDD;
  return vel_cov_mtx;
}

void init_controller_and_estimator_constants() {
  constantsASTRA.g = 9.8015;
  constantsASTRA.m = 1.1690;
  constantsASTRA.mag << 0.38535, 0.03198, -0.9222;
  constantsASTRA.rTB = 0.26f; // TODO: Check if true
  constantsASTRA.MaxThrust = 16.0f;
  constantsASTRA.J << 0.067f, 0.067f, 0.02f;
  constantsASTRA.K_I = Matrix3_3::Identity();
  constantsASTRA.K_I.diagonal() << 0.15f, 0.15f, 0.05f;
  constantsASTRA.Q = Q_gen();

  // Gain Matrices (TODO: Check values)
  constantsASTRA.K_trans << 
    9.982694e-01, 0, 0, 1.730052e+00, 0, 0, 
    0, 9.982694e-01, 0, 0, 1.730052e+00, 0, 
    0, 0, 2.988404e+00, 0, 0, 3.861005e+00;  
  constantsASTRA.K_rot << 
    -5.867653e-01, -7.171442e-05, -6.699845e-06, -2.179070e-01, 1.061392e-06, 5.317069e-06, 
    7.171184e-05, -5.867654e-01, -1.629788e-05, -1.043975e-06, -2.179071e-01, 2.393766e-06, 
    -1.919781e-05, -6.319698e-06, 8.237270e-01, -4.350950e-06, -1.350637e-06, 2.299779e-01; 

  constantsASTRA.L_Thrust << 5.616296e-03, 0, 0, 2.005830e-03, 0, 0, 0, 4.727018e-03, 4.929804e-06, 0, 2.010921e-03, 2.808143e-08, 0, 5.023681e-06, 5.088132e-03, 0, -1.601761e-07, 2.006986e-03,
      1.889238e-04, -9.886904e-10, 0, 1.123218e-02, 0, 0, 0, 9.157307e-04, -2.406860e-05, 0, 9.845408e-03, -4.535823e-07, 0, 3.624535e-06, 2.974015e-04, 0, -3.739723e-07, 9.309009e-03, 5.350414e-04,
      -2.769941e-09, 0, 1.575495e-02, 1.030446e-10, 0, 0, 2.273589e-03, -5.977943e-05, 0, 1.205917e-02, -1.064738e-06, 0, 8.297828e-06, 6.894103e-04, 0, -8.723946e-07, 1.079742e-02;
  constantsASTRA.L_Att << 6.808373e-03, 0, 0, 0, 8.123445e-03, 0, 0, 0, 8.123349e-03, 5.788234e-03, 0, 0, 0, 8.242795e-03, 0, 0, 0, 8.194591e-03;

  constantsASTRA.R = Matrix6_6::Zero();
  constantsASTRA.R.block<3, 3>(0, 0) = Matrix3_3::Identity() * 0.002f;
  constantsASTRA.R.block<3, 3>(3, 3) = Matrix3_3::Identity() * 0.050f;

  P = Matrix15_15::Identity();
  P.block<3, 3>(0, 0) = Matrix3_3::Identity() * 0.05f;
  P.block<3, 3>(3, 3) = Matrix3_3::Identity() * 1.0f;
  P.block<3, 3>(6, 6) = Matrix3_3::Identity() * 0.1f;
  P.block<3, 3>(9, 9) = Matrix3_3::Identity() * 0.001f;
  P(12, 12) = 0.0004f;
  P(13, 13) = 0.0004f;
  P(14, 14) = 0.01f;

  x_est = Vector19::Zero();
  x_est[0] = 1;
  init_flight_covariance();
  lastEMA = Vector2::Zero();
  dnf_X = Matrix9_4::Ones();
  dnf_Y = Matrix9_4::Ones();
  xhat_pos.setZero();
  xhat_att.setZero();
  last_U_cmd.setZero();
  last_thrust = 0;
  last_GND = true;

  ASTRAv2_Controller_reset();
}

Controller_Output get_controller_output(Controller_Input ci, const Vector15 &X_trg, const Vector4 &U_ff, float ideal_dT, float loop_dT, Controller_Internals *cs) {
  Vector15 z;
  // clang-format off
  z << ci.imu.accel_x, ci.imu.accel_y, ci.imu.accel_z,
       ci.imu.gyro_yaw, ci.imu.gyro_pitch, ci.imu.gyro_roll,
       ci.mag.mag_x, ci.mag.mag_y, ci.mag.mag_z,
       ci.gps.pos.north, ci.gps.pos.west, ci.gps.pos.up,
       ci.gps.vel.north, ci.gps.vel.west, ci.gps.vel.up;
  // clang-format on

  Vector9 imu = z.segment<9>(0);
  if (last_thrust < 1) { // prevent div by 0 in filter
    last_thrust = 9.8;
  }
  Vector9 filt_imu = DigitalNF(imu, ci.GND_val, last_thrust / constantsASTRA.MaxThrust * 100.0, ideal_dT, dnf_X, dnf_Y);
  z.segment<9>(0) = filt_imu;

  for (int i = 0; i < 9; i++) {
    cs->filter_out[i] = z(i);
  }

  if (!ci.GND_val && last_GND) { // we left GND this frame
    init_flight_covariance();
    ASTRAv2_Controller_reset(); // reset integral gains in the controller itself
  }

  Matrix3_3 gps_vel_covar = get_vel_cov_mtx(ci.gps);
  Matrix3_3 gps_pos_covar = get_pos_cov_mtx(ci.gps);

  if (ci.GND_val) {
    x_est = GroundEstimator(x_est, constantsASTRA, z, loop_dT, P, ci.new_mag_packet, ci.new_gps_packet, gps_vel_covar, gps_pos_covar);
  } else {
    x_est = FlightEstimator2(x_est, constantsASTRA, z, loop_dT, Flight_P, ci.new_mag_packet, ci.new_gps_packet, gps_vel_covar, gps_pos_covar);
  }

  // CHANGE: Bias corrected rates
  Vector3 omega_body = z.segment<3>(3) - x_est.segment<3>(10);
  Vector15 X_ctrl;
  X_ctrl << x_est.segment<10>(0), omega_body, 0.0f, 0.0f;

  // Run disturbance observers
  Vector3 a_dist = LESO_Position(ci.GND_val, X_ctrl, last_U_cmd, constantsASTRA.L_Thrust, constantsASTRA, loop_dT, xhat_pos);
  Vector3 u_dist = LESO_Attitude(ci.GND_val, X_ctrl, last_U_cmd, constantsASTRA.L_Att, constantsASTRA, loop_dT, xhat_att);
  Vector6 Dist_LESO;
  Dist_LESO << a_dist, u_dist;
  last_X_ctrl = X_ctrl; // abort harness: copies only, no effect on the estimator or controller
  last_a_dist = a_dist;

  // Main control loop
  Vector4 raw_co = ASTRAv2_Controller(ci.GND_val, X_ctrl, X_trg, U_ff, constantsASTRA, Dist_LESO, loop_dT, cs);
  last_U_cmd = raw_co;
  if (ci.GND_val) {
    raw_co = Vector4::Zero();
  }
  last_thrust = raw_co(2);
  last_GND = ci.GND_val;

  float co_gimbal_yaw = raw_co(0) * 180 / M_PI;
  float co_gimbal_pitch = raw_co(1) * 180 / M_PI;

  float ALPHA = ideal_dT / tau;
  float new_gimbal_yaw = ALPHA * co_gimbal_yaw + (1.0f - ALPHA) * lastEMA(0);
  float new_gimbal_pitch = ALPHA * co_gimbal_pitch + (1.0f - ALPHA) * lastEMA(1);

  lastEMA(0) = new_gimbal_yaw;
  lastEMA(1) = new_gimbal_pitch;

  Controller_Output co;
  co.gimbal_yaw_deg = new_gimbal_yaw;
  co.gimbal_pitch_deg = new_gimbal_pitch;
  co.thrust_N = raw_co(2);
  co.roll_rad_sec_squared = raw_co(3);
  co.dist_accel_x = a_dist(0);
  co.dist_accel_y = a_dist(1);
  co.dist_accel_z = a_dist(2);
  co.dist_rate_x  = u_dist(0);
  co.dist_rate_y  = u_dist(1);
  co.dist_rate_z  = u_dist(2);

  // Add the full estimated state
  cs->x_est.q_vec_w = x_est(0);
  cs->x_est.q_vec_x = x_est(1);
  cs->x_est.q_vec_y = x_est(2);
  cs->x_est.q_vec_z = x_est(3);
  cs->x_est.est_pos_north = x_est(4);
  cs->x_est.est_pos_west = x_est(5);
  cs->x_est.est_pos_up = x_est(6);
  cs->x_est.est_vel_north = x_est(7);
  cs->x_est.est_vel_west = x_est(8);
  cs->x_est.est_vel_up = x_est(9);
  cs->x_est.gyro_bias_yaw = x_est(10);
  cs->x_est.gyro_bias_pitch = x_est(11);
  cs->x_est.gyro_bias_roll = x_est(12);
  cs->x_est.accel_bias_x = x_est(13);
  cs->x_est.accel_bias_y = x_est(14);
  cs->x_est.accel_bias_z = x_est(15);
  cs->x_est.mag_bias_x = x_est(16);
  cs->x_est.mag_bias_y = x_est(17);
  cs->x_est.mag_bias_z = x_est(18);
  return co;
}

const Vector15 &get_controller_x_est() { return last_X_ctrl; }

const Vector3 &get_position_disturbance() { return last_a_dist; }
} // namespace ControllerAndEstimator
