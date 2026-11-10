// Rotational Disturbance Observer
// Pablo Plata 10/13/26

#include "matlab_funcs.h"

Vector3 LESO_Attitude(bool GND, const Vector15 &X_est, const Vector4 &U_cmd, const Matrix6_3 &L_Att, const constantsASTRA_t &constants, float dT, Vector6 &xhat) {

    // the following vars are persistent: xhat
    if (dT <= 0.0f || GND == 1) {
        xhat.setZero();
        return Vector3::Zero();
    }

    // Setup Dynamics & State Propagation
    Vector3 y_meas = X_est.segment<3>(10);

    // Nominal thrust force vector in body frame
    float theta = U_cmd[0];
    float phi = U_cmd[1];
    float thrust = U_cmd[2];
    float roll = U_cmd[3];

    Vector3 ThrustVec_B = thrust * (Vector3() << cos(theta) * sin(phi), -sin(theta), cos(theta) * cos(phi)).finished();

    // Commanded Torque via RBD
    Vector3 leverArm;
    leverArm << 0.0f, 0.0f, -constants.rTB;
    Vector3 tau_cmd = leverArm.cross(ThrustVec_B) + Vector3(0.0f, 0.0f, roll);
    Vector3 angAccel = tau_cmd.cwiseQuotient(constants.J);

    // Gyroscopic Coupling
    Vector3 J_omega = y_meas.cwiseProduct(constants.J);
    Vector3 gyroCross = y_meas.cross(J_omega);
    Vector3 gyroCouple = gyroCross.cwiseQuotient(constants.J);

    // Discrete STMs
    Matrix6_6 A_LESO_Att = Matrix6_6::Identity();
    A_LESO_Att.block<3, 3>(0, 3) = Matrix3_3::Identity() * dT;

    Matrix6_3 B_LESO_Att = Matrix6_3::Zero();
    B_LESO_Att.block<3, 3>(0, 0) = Matrix3_3::Identity() * dT;

    // State Prediction
    Vector6 xhat_pred = A_LESO_Att * xhat + B_LESO_Att * (angAccel - gyroCouple);

    // State Update
    Vector3 y_hat = xhat_pred.segment<3>(0);
    xhat = xhat_pred + L_Att * (y_meas - y_hat);

    // Disturbance Acceleration
    Vector3 a_dist = xhat.segment<3>(3);
    return a_dist;
}