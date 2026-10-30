#define S_FUNCTION_NAME  ASTRAv2_SMEX
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"

// Flight software headers
#include "controller_and_estimator.h"
#include "matlab_funcs.h"

static void mdlInitializeSizes(SimStruct *S) {
    // Zero parameters needed for MVP (assumes 1 kHz / 0.001 s loop)
    ssSetNumSFcnParams(S, 0);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; // Parameter mismatch
    }

    // Input Ports
    // Port 0: Sensor Measurements [accel(3); gyro(3); mag(3); gps_pos(3); gps_vel(3)] = 15 floats
    // Port 1: Guidance Target State X_trg = 15 floats
    // Port 2: Feedforward Command U_ff = 4 floats
    // Port 3: Status Flags [GND_val; new_mag_packet; new_gps_packet] = 3 booleans
    if (!ssSetNumInputPorts(S, 4)) return;

    // Port 0: Sensors
    ssSetInputPortWidth(S, 0, 15);
    ssSetInputPortDataType(S, 0, SS_SINGLE);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortRequiredContiguous(S, 0, 1);

    // Port 1: Target State
    ssSetInputPortWidth(S, 1, 15);
    ssSetInputPortDataType(S, 1, SS_SINGLE);
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    ssSetInputPortRequiredContiguous(S, 1, 1);

    // Port 2: Feedforward
    ssSetInputPortWidth(S, 2, 4);
    ssSetInputPortDataType(S, 2, SS_SINGLE);
    ssSetInputPortDirectFeedThrough(S, 2, 1);
    ssSetInputPortRequiredContiguous(S, 2, 1);

    // Port 3: Flags
    ssSetInputPortWidth(S, 3, 3);
    ssSetInputPortDataType(S, 3, SS_BOOLEAN);
    ssSetInputPortDirectFeedThrough(S, 3, 1);
    ssSetInputPortRequiredContiguous(S, 3, 1);

    // Output Ports
    // Port 0: Actuator Commands U_cmd [pitch_deg; yaw_deg; thrust_N; roll_accel] = 4 floats
    // Port 1: Estimated State X_est [q(4); pos(3); vel(3); biases(9)] = 19 floats
    if (!ssSetNumOutputPorts(S, 2)) return;

    // Port 0: Commands
    ssSetOutputPortWidth(S, 0, 4);
    ssSetOutputPortDataType(S, 0, SS_SINGLE);

    // Port 1: Telemetry / State Estimate
    ssSetOutputPortWidth(S, 1, 19);
    ssSetOutputPortDataType(S, 1, SS_SINGLE);

    // Timing & Engine Options
    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    // Support Simulink Accelerator
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE |
                    SS_OPTION_USE_TLC_WITH_ACCELERATOR);
}

static void mdlInitializeSampleTimes(SimStruct *S) {
    // 1 kHz discrete sample time (0.001 s)
    ssSetSampleTime(S, 0, 0.001);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S);
}

#define MDL_START
#if defined(MDL_START)
static void mdlStart(SimStruct *S) {
    // Initialize controller and estimator states and constants
    ControllerAndEstimator::init_controller_and_estimator_constants();
}
#endif

static void mdlOutputs(SimStruct *S, int_T tid) {
    // Acquire input pointers
    const real32_T  *y_meas = (const real32_T*)  ssGetInputPortSignal(S, 0);
    const real32_T  *x_trg  = (const real32_T*)  ssGetInputPortSignal(S, 1);
    const real32_T  *u_ff   = (const real32_T*)  ssGetInputPortSignal(S, 2);
    const boolean_T *flags  = (const boolean_T*) ssGetInputPortSignal(S, 3);

    // Acquire output pointers (using ssGetOutputPortSignal for typed single-precision)
    real32_T *u_cmd = (real32_T*) ssGetOutputPortSignal(S, 0);
    real32_T *x_est = (real32_T*) ssGetOutputPortSignal(S, 1);

    // Build Controller_Input struct
    Controller_Input ci{};
    ci.GND_val        = (flags[0] != 0);
    ci.new_mag_packet = (flags[1] != 0);
    ci.new_gps_packet = (flags[2] != 0);

    // IMU measurements
    ci.imu.accel_x    = y_meas[0];
    ci.imu.accel_y    = y_meas[1];
    ci.imu.accel_z    = y_meas[2];
    ci.imu.gyro_yaw   = y_meas[3];
    ci.imu.gyro_pitch = y_meas[4];
    ci.imu.gyro_roll  = y_meas[5];

    // Magnetometer measurements
    ci.mag.mag_x = y_meas[6];
    ci.mag.mag_y = y_meas[7];
    ci.mag.mag_z = y_meas[8];

    // GPS measurements
    ci.gps.pos.north = y_meas[9];
    ci.gps.pos.west  = y_meas[10];
    ci.gps.pos.up    = y_meas[11];
    ci.gps.vel.north = y_meas[12];
    ci.gps.vel.west  = y_meas[13];
    ci.gps.vel.up    = y_meas[14];

    // Isotropic GPS covariances
    ci.gps.posCovNN = 1.0f; ci.gps.posCovNE = 0.0f; ci.gps.posCovND = 0.0f;
    ci.gps.posCovEE = 1.0f; ci.gps.posCovED = 0.0f; ci.gps.posCovDD = 1.0f;
    ci.gps.velCovNN = 1.0f; ci.gps.velCovNE = 0.0f; ci.gps.velCovND = 0.0f;
    ci.gps.velCovEE = 1.0f; ci.gps.velCovED = 0.0f; ci.gps.velCovDD = 1.0f;

    // Zero-copy wrapping of input vectors into Eigen matrices
    Eigen::Map<const Vector15> X_trg_vec(x_trg);
    Eigen::Map<const Vector4>  U_ff_vec(u_ff);

    // Step the C++ flight controller and estimator
    const float dt = 0.001f;
    Controller_Internals cs{};
    Controller_Output co = ControllerAndEstimator::get_controller_output(ci, X_trg_vec, U_ff_vec, dt, dt, &cs);

    // Assign commands to Output Port 0
    u_cmd[0] = co.gimbal_pitch_deg;
    u_cmd[1] = co.gimbal_yaw_deg;
    u_cmd[2] = co.thrust_N;
    u_cmd[3] = co.roll_rad_sec_squared;

    // Assign estimated state telemetry to Output Port 1
    x_est[0]  = cs.x_est.q_vec_w;
    x_est[1]  = cs.x_est.q_vec_x;
    x_est[2]  = cs.x_est.q_vec_y;
    x_est[3]  = cs.x_est.q_vec_z;
    x_est[4]  = cs.x_est.est_pos_north;
    x_est[5]  = cs.x_est.est_pos_west;
    x_est[6]  = cs.x_est.est_pos_up;
    x_est[7]  = cs.x_est.est_vel_north;
    x_est[8]  = cs.x_est.est_vel_west;
    x_est[9]  = cs.x_est.est_vel_up;
    x_est[10] = cs.x_est.gyro_bias_yaw;
    x_est[11] = cs.x_est.gyro_bias_pitch;
    x_est[12] = cs.x_est.gyro_bias_roll;
    x_est[13] = cs.x_est.accel_bias_x;
    x_est[14] = cs.x_est.accel_bias_y;
    x_est[15] = cs.x_est.accel_bias_z;
    x_est[16] = cs.x_est.mag_bias_x;
    x_est[17] = cs.x_est.mag_bias_y;
    x_est[18] = cs.x_est.mag_bias_z;
}

static void mdlTerminate(SimStruct *S) {
    // Teardown (no-op)
}

#ifdef MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
