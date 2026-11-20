#pragma once
#ifndef ABORT_CODE_H
#define ABORT_CODE_H

#include "matlab_funcs.h"

// Guidance Constants
constexpr float ABORT_HOVER_ALT = 1.0f;       
constexpr float ABORT_R_MAX = 3.0f;           
constexpr float ABORT_DT_SAMPLE = 0.2f;       
constexpr int ABORT_MAX_POINTS = 100;         
constexpr float ABORT_MIN_SEG_TIME = 0.2f;    
constexpr float ABORT_MAX_T1 = 15.0f;         
constexpr int ABORT_MAX_DILATE_ITERS = 10;    // Iterations for time allocation
constexpr int ABORT_EVAL_GRID = 41;           

// Vehicle Limits
constexpr float ABORT_THRUST_MARGIN = 0.10f;
constexpr float ABORT_GIMBAL_MARGIN = 0.10f;
constexpr float ABORT_MIN_THROTTLE = 0.25f;
constexpr float ABORT_DEG2RAD = 3.1416f / 180.0f;
constexpr float ABORT_MAX_GIMBAL_RAD = 10.0f * ABORT_DEG2RAD;
constexpr float ABORT_MAX_ROLL = 5.0f;
constexpr float ABORT_MAX_TILT_RAD = 30.0f * ABORT_DEG2RAD;
constexpr float ABORT_MAX_VEL_XY = 4.0f;
constexpr float ABORT_MAX_VEL_Z = 5.0f;
constexpr float ABORT_LAND_VEL_Z = 1.0f;
constexpr float ABORT_DEF_ACCEL_XY = 2.5f;
constexpr float ABORT_DEF_ACCEL_Z = 3.5f;

struct AbortWrenchPoint {
    Vector3 F_body;
    Vector3 Tau_body;
};

struct AbortTrajectory {
    int num_points;
    float time[ABORT_MAX_POINTS];
    float states[ABORT_MAX_POINTS][15];
    float inputs[ABORT_MAX_POINTS][4];
};

void flatnessInvertPoint(const Vector3 &P_pos, const Vector3 &P_vel, const Vector3 &P_acc, const Vector3 &P_jerk, const Vector3 &P_snap, float rho, const constantsASTRA_t &constants,
                         const Vector4 *q_prev, Vector15 &X_out, AbortWrenchPoint &Wrench_out, Vector3 &AccelCg_out);

Vector4 allocateControlsPoint(const AbortWrenchPoint &Wrench, float rTB);

void computeAbort(float t_now, const Vector15 &X_est, float last_thrust, const constantsASTRA_t &constants, AbortTrajectory &traj);

#endif // ABORT_CODE_H
