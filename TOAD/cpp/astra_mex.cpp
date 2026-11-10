#include "mex.h"
#include <cstring>
#include "controller_and_estimator.h"

using namespace ControllerAndEstimator;

//helpers
static const mxArray *sub(const mxArray *s, const char *name) {
  if (!mxIsStruct(s)) mexErrMsgIdAndTxt("astra:type", "Expected a struct while reading '%s'", name);
  const mxArray *f = mxGetField(s, 0, name);
  if (!f) mexErrMsgIdAndTxt("astra:field", "Missing struct field '%s'", name);
  return f;
}
static float num(const mxArray *s, const char *name) {
  const mxArray *f = sub(s, name);
  if (!(mxIsNumeric(f) || mxIsLogical(f)) || mxGetNumberOfElements(f) != 1)
    mexErrMsgIdAndTxt("astra:type", "Field '%s' must be a numeric/logical scalar", name);
  return (float)mxGetScalar(f);
}
static float numOpt(const mxArray *s, const char *name, float def) {
  const mxArray *f = mxGetField(s, 0, name);
  return f ? (float)mxGetScalar(f) : def;
}
static void put(mxArray *s, const char *name, mxArray *val) {
  mxAddField(s, name);
  mxSetField(s, 0, name, val);
}
static mxArray *vec(const float *v, int n) {
  mxArray *a = mxCreateDoubleMatrix(n, 1, mxREAL);
  double *p = mxGetPr(a);
  for (int i = 0; i < n; i++) p[i] = v[i];
  return a;
}
static mxArray *scalar(double v) { return mxCreateDoubleScalar(v); }

//struct
static Controller_Input readInput(const mxArray *s) {
  Controller_Input ci{};
  ci.GND_val = num(s, "GND_val") != 0;
  ci.new_mag_packet = num(s, "new_mag_packet") != 0;
  ci.new_gps_packet = num(s, "new_gps_packet") != 0;

  const mxArray *imu = sub(s, "imu");
  ci.imu.accel_x = num(imu, "accel_x");
  ci.imu.accel_y = num(imu, "accel_y");
  ci.imu.accel_z = num(imu, "accel_z");
  ci.imu.gyro_yaw = num(imu, "gyro_yaw");
  ci.imu.gyro_pitch = num(imu, "gyro_pitch");
  ci.imu.gyro_roll = num(imu, "gyro_roll");

  const mxArray *mag = sub(s, "mag");
  ci.mag.mag_x = num(mag, "mag_x");
  ci.mag.mag_y = num(mag, "mag_y");
  ci.mag.mag_z = num(mag, "mag_z");

  const mxArray *gps = sub(s, "gps");
  const mxArray *pos = sub(gps, "pos");
  const mxArray *vel = sub(gps, "vel");
  ci.gps.pos.north = num(pos, "north");
  ci.gps.pos.west = num(pos, "west");
  ci.gps.pos.up = num(pos, "up");
  ci.gps.vel.north = num(vel, "north");
  ci.gps.vel.west = num(vel, "west");
  ci.gps.vel.up = num(vel, "up");

  ci.gps.posCovNN = num(gps, "posCovNN");
  ci.gps.posCovNE = num(gps, "posCovNE");
  ci.gps.posCovND = num(gps, "posCovND");
  ci.gps.posCovEE = num(gps, "posCovEE");
  ci.gps.posCovED = num(gps, "posCovED");
  ci.gps.posCovDD = num(gps, "posCovDD");
  ci.gps.velCovNN = num(gps, "velCovNN");
  ci.gps.velCovNE = num(gps, "velCovNE");
  ci.gps.velCovND = num(gps, "velCovND");
  ci.gps.velCovEE = num(gps, "velCovEE");
  ci.gps.velCovED = num(gps, "velCovED");
  ci.gps.velCovDD = num(gps, "velCovDD");

  ci.target_pos_north = numOpt(s, "target_pos_north", 0);
  ci.target_pos_west = numOpt(s, "target_pos_west", 0);
  ci.target_pos_up = numOpt(s, "target_pos_up", 0);
  return ci;
}

static void readVec(const mxArray *a, const char *name, int n, float *out) {
  if (!mxIsDouble(a) || mxIsComplex(a) || (int)mxGetNumberOfElements(a) != n)
    mexErrMsgIdAndTxt("astra:size", "%s must be a real double vector with %d elements", name, n);
  const double *p = mxGetPr(a);
  for (int i = 0; i < n; i++) out[i] = (float)p[i];
}

//actual MEX stuff
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs < 1 || !mxIsChar(prhs[0]))
    mexErrMsgIdAndTxt("astra:usage", "First argument must be a command: 'init' or 'step'");

  char cmd[16];
  mxGetString(prhs[0], cmd, sizeof(cmd));

  if (std::strcmp(cmd, "init") == 0) {
    init_controller_and_estimator_constants();
    return;
  }

  if (std::strcmp(cmd, "step") == 0) {
    if (nrhs != 6)
      mexErrMsgIdAndTxt("astra:usage", "Usage: [co,int] = astra_mex('step', ci, X_trg, U_ff, ideal_dT, loop_dT)");

    Controller_Input ci = readInput(prhs[1]);

    float xt[15], uf[4];
    readVec(prhs[2], "X_trg", 15, xt);
    readVec(prhs[3], "U_ff", 4, uf);
    Vector15 X_trg;
    Vector4 U_ff;
    for (int i = 0; i < 15; i++) X_trg(i) = xt[i];
    for (int i = 0; i < 4; i++) U_ff(i) = uf[i];

    float ideal_dT = (float)mxGetScalar(prhs[4]);
    float loop_dT = (float)mxGetScalar(prhs[5]);

    Controller_Internals cs{};
    Controller_Output co = get_controller_output(ci, X_trg, U_ff, ideal_dT, loop_dT, &cs);

    //controller output
    plhs[0] = mxCreateStructMatrix(1, 1, 0, nullptr);
    put(plhs[0], "thrust_N", scalar(co.thrust_N));
    put(plhs[0], "roll_rad_sec_squared", scalar(co.roll_rad_sec_squared));
    put(plhs[0], "gimbal_pitch_deg", scalar(co.gimbal_pitch_deg));
    put(plhs[0], "gimbal_yaw_deg", scalar(co.gimbal_yaw_deg));

    //controller
    if (nlhs > 1) {
      plhs[1] = mxCreateStructMatrix(1, 1, 0, nullptr);

      //x hat
      mxArray *xe = mxCreateStructMatrix(1, 1, 0, nullptr);
      const Controller_X_Est &x = cs.x_est;
      put(xe, "q_vec_w", scalar(x.q_vec_w));
      put(xe, "q_vec_x", scalar(x.q_vec_x));
      put(xe, "q_vec_y", scalar(x.q_vec_y));
      put(xe, "q_vec_z", scalar(x.q_vec_z));
      put(xe, "est_pos_north", scalar(x.est_pos_north));
      put(xe, "est_pos_west", scalar(x.est_pos_west));
      put(xe, "est_pos_up", scalar(x.est_pos_up));
      put(xe, "est_vel_north", scalar(x.est_vel_north));
      put(xe, "est_vel_west", scalar(x.est_vel_west));
      put(xe, "est_vel_up", scalar(x.est_vel_up));
      put(xe, "gyro_bias_yaw", scalar(x.gyro_bias_yaw));
      put(xe, "gyro_bias_pitch", scalar(x.gyro_bias_pitch));
      put(xe, "gyro_bias_roll", scalar(x.gyro_bias_roll));
      put(xe, "accel_bias_x", scalar(x.accel_bias_x));
      put(xe, "accel_bias_y", scalar(x.accel_bias_y));
      put(xe, "accel_bias_z", scalar(x.accel_bias_z));
      put(xe, "mag_bias_x", scalar(x.mag_bias_x));
      put(xe, "mag_bias_y", scalar(x.mag_bias_y));
      put(xe, "mag_bias_z", scalar(x.mag_bias_z));
      put(plhs[1], "x_est", xe);

      put(plhs[1], "filter_out", vec(cs.filter_out, 9));
      put(plhs[1], "VelTarget", vec(cs.VelTarget, 3));
      put(plhs[1], "AccelTarget", vec(cs.AccelTarget, 3));
      put(plhs[1], "TargetAtt", vec(cs.TargetAtt, 4));
      put(plhs[1], "VelErrorI", vec(cs.VelErrorI, 3));
      put(plhs[1], "AttErrorI", vec(cs.AttErrorI, 3));
    }
    return;
  }

  mexErrMsgIdAndTxt("astra:usage", "Unknown command '%s'", cmd);
}