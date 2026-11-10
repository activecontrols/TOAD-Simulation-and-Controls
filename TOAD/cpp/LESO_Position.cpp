// Translational Disturbance Observer
// Pablo Plata 10/13/26

#include "matlab_funcs.h"

Vector3 LESO_Position(bool GND, const Vector15 &X_est, const Vector4 &U_cmd, const Matrix9_6 &L_Thrust, const constantsASTRA_t &constants, float dT, Vector9 &xhat) {

    // the following vars are persistent: xhat
    if (dT <= 0.0f || GND == 1) {
        xhat.setZero();
        return Vector3::Zero();
    }

    // Setup Dynamics & State Propagation
    Vector6 y_meas;
    y_meas << X_est.segment<3>(4), X_est.segment<3>(7);

    float Mass = constants.m + X_est[13] + X_est[14];
    Matrix3_3 C_BI = quatRot(X_est.segment<4>(0));
    Matrix3_3 C_IB = C_BI.transpose();

    // Nominal thrust force vector in body frame
    float theta = U_cmd[0];
    float phi = U_cmd[1];
    float thrust = U_cmd[2];
    Vector3 ThrustVec_B = thrust * (Vector3() << cos(theta) * sin(phi), -sin(theta), cos(theta) * cos(phi)).finished();

    // Nominal acceleration in inertial frame (inc. gravity correction)
    Vector3 g_vec = (Vector3() << 0, 0, -constants.g).finished();
    Vector3 a_nominal_I = (C_IB * ThrustVec_B) / Mass + g_vec;

    // Discrete STMs
    Matrix9_9 A_LESO_Thr = Matrix9_9::Identity();
    A_LESO_Thr.block<3, 3>(0, 3) = Matrix3_3::Identity() * dT;
    A_LESO_Thr.block<3, 3>(0, 6) = Matrix3_3::Identity() * (0.5f * dT * dT);
    A_LESO_Thr.block<3, 3>(3, 6) = Matrix3_3::Identity() * dT;

    Matrix9_3 B_LESO_Thr = Matrix9_3::Zero();
    B_LESO_Thr.block<3, 3>(0, 0) = Matrix3_3::Identity() * (0.5f * dT * dT);
    B_LESO_Thr.block<3, 3>(3, 0) = Matrix3_3::Identity() * dT;

    // State Prediction
    Vector9 xhat_pred = A_LESO_Thr * xhat + B_LESO_Thr * a_nominal_I;

    // State Update
    Vector6 y_hat = xhat_pred.segment<6>(0);
    xhat = xhat_pred + L_Thrust * (y_meas - y_hat);

    // Disturbance Acceleration
    Vector3 a_dist = xhat.segment<3>(6);
    return a_dist;
}