#pragma once
#ifndef MATLAB_FUNCS_H
#define MATLAB_FUNCS_H

#include "astra_structs.h"
#include <ArduinoEigenDense.h>

// Note - using matlab funcs based off of ASTRA simulation: 3a2af9e756f0981da411d93f59ee6b421d5fa234
// TODO - update matlab func commit

using Vector2 = Eigen::Matrix<float, 2, 1>;
using Vector3 = Eigen::Matrix<float, 3, 1>;
using Vector4 = Eigen::Matrix<float, 4, 1>;
using Vector6 = Eigen::Matrix<float, 6, 1>;
using Vector9 = Eigen::Matrix<float, 9, 1>;
using Vector12 = Eigen::Matrix<float, 12, 1>;
using Vector13 = Eigen::Matrix<float, 13, 1>;
using Vector15 = Eigen::Matrix<float, 15, 1>;
using Vector16 = Eigen::Matrix<float, 16, 1>;
using Vector18 = Eigen::Matrix<float, 18, 1>;
using Vector19 = Eigen::Matrix<float, 19, 1>;
using Matrix2_2 = Eigen::Matrix<float, 2, 2>;
using Matrix3_3 = Eigen::Matrix<float, 3, 3>;
using Matrix3_6 = Eigen::Matrix<float, 3, 6>;
using Matrix6_6 = Eigen::Matrix<float, 6, 6>;
using Matrix6_3 = Eigen::Matrix<float, 6, 3>;
using Matrix3_9 = Eigen::Matrix<float, 3, 9>;
using Matrix3_12 = Eigen::Matrix<float, 3, 12>;
using Matrix4_4 = Eigen::Matrix<float, 4, 4>;
using Matrix6_9 = Eigen::Matrix<float, 6, 9>;
using Matrix6_12 = Eigen::Matrix<float, 6, 12>;
using Matrix9_3 = Eigen::Matrix<float, 9, 3>;
using Matrix9_6 = Eigen::Matrix<float, 9, 6>;
using Matrix9_9 = Eigen::Matrix<float, 9, 9>;
using Matrix12_3 = Eigen::Matrix<float, 12, 3>;
using Matrix12_6 = Eigen::Matrix<float, 12, 6>;
using Matrix12_12 = Eigen::Matrix<float, 12, 12>;
using Matrix6_18 = Eigen::Matrix<float, 6, 18>;
using Matrix18_6 = Eigen::Matrix<float, 18, 6>;
using Matrix18_18 = Eigen::Matrix<float, 18, 18>;

typedef struct {
    // Kalman Filt
    Vector3 mag;
    Matrix18_18 Q;
    Matrix6_6 R;

    // Controls
    Matrix3_6 K_trans;
    Matrix3_6 K_rot;
    Matrix3_3 K_I;

    // LESOs
    Matrix9_6 L_Thrust;
    Matrix6_3 L_Att;

    // Vehicle Params
    float g;
    float m;
    float rTB;
    float MaxThrust;
    Vector3 J;

} constantsASTRA_t;

Matrix3_3 zetaCross(Vector3 zeta);
Matrix3_3 quatRot(Vector4 q);
Matrix4_4 HamiltonianProd(Vector4 q);
Vector19 GroundEstimator(Vector19 x_est, constantsASTRA_t constantsASTRA, Vector15 z, float dT, Matrix18_18 &P, bool new_mag_packet, bool new_gps_packet, Matrix3_3 gps_vel_covar,
                         Matrix3_3 gps_pos_covar);
Vector19 FlightEstimator(Vector19 x_est, constantsASTRA_t constantsASTRA, Vector15 z, float dT, Matrix9_9 &P, bool new_gps_packet, Matrix3_3 gps_vel_covar, Matrix3_3 gps_pos_covar);
Vector19 FlightEstimator2(Vector19 x_est, constantsASTRA_t constantsASTRA, Vector15 z, float dT, Matrix12_12 &P, bool new_mag_packet, bool new_gps_packet, Matrix3_3 gps_vel_covar, Matrix3_3 gps_pos_covar);
Vector16 StateAUG(Vector13 XKF, Vector3 G);
Vector4 DCM_Quat_Conversion(Matrix3_3 R);
Vector9 DigitalNF(Vector9 IN, bool GND, float THRUST_PERC, float dT, Matrix9_4 &X, Matrix9_4 &Y);
Matrix18_18 Q_gen();
void ASTRAv2_Controller_reset();
Vector4 ASTRAv2_Controller(bool GND, const Vector15 &X_est, const Vector15 &X_trg, const Vector4 &U_ff, const constantsASTRA_t &constantsTOAD, const Vector6 &Dist_LESO, float dT,
                           Controller_Internals *cs);
Vector3 LESO_Position(bool GND, const Vector15 &X_est, const Vector4 &U_cmd, const Matrix9_6 &L_Thrust, const constantsASTRA_t &constants, float dT, Vector9 &xhat);
Vector3 LESO_Attitude(bool GND, const Vector15 &X_est, const Vector4 &U_cmd, const Matrix6_3 &L_Att, const constantsASTRA_t &constants, float dT, Vector6 &xhat);

// call like so: matrixExpPade6<Matrix9_9>
template <typename MatrixN_N> MatrixN_N matrixExpPade6(MatrixN_N A) {
    using Scalar = typename MatrixN_N::Scalar;
    const Scalar b[] = {1.0, 0.5, 0.12, 0.0183333333333, 0.00199275362319, 0.000160590438368, 0.00000939085239292};

    MatrixN_N A2 = A * A;
    MatrixN_N A4 = A2 * A2;
    MatrixN_N A6 = A4 * A2;

    MatrixN_N U = A * (b[1] * MatrixN_N::Identity(A.rows(), A.cols()) + b[3] * A2 + b[5] * A4);
    U += A6 * b[6];

    MatrixN_N V = b[0] * MatrixN_N::Identity(A.rows(), A.cols()) + b[2] * A2 + b[4] * A4 + b[6] * A6;

    MatrixN_N numer = V + U;
    MatrixN_N denom = V - U;

    return denom.inverse() * numer;
}

#endif // MATLAB_FUNCS_H
