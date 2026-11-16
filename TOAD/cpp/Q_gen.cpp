#include "matlab_funcs.h"

Matrix18_18 Q_gen(float dt) {
    if (dt <= 0.0f) {
        return Matrix18_18::Zero();
    }

    const float gyro_cov = 5e-6f;
    const float gyro_bias_cov = 1e-6f;
    const float accel_proc_cov = 1e-2f;
    const float accel_bias_cov = 1e-8f;
    // Ornstein-Uhlenbeck mag bias slack: tau_m = 60 s, sigma_ss = 0.015 raw units -> q = 2*sigma_ss^2/tau_m
    const float mag_bias_cov = 7.5e-6f;

    const float dt2 = dt * dt;
    const float dt3 = dt2 * dt;
    const float dt4 = dt3 * dt;
    const float dt5 = dt4 * dt;

    Matrix3_3 gyro_cov_mat = gyro_cov * Matrix3_3::Identity();
    Matrix3_3 gyro_bias_cov_mat = gyro_bias_cov * Matrix3_3::Identity();
    Matrix3_3 accel_cov_mat = accel_proc_cov * Matrix3_3::Identity();
    Matrix3_3 accel_bias_cov_mat = accel_bias_cov * Matrix3_3::Identity();
    Matrix3_3 mag_bias_cov_mat = mag_bias_cov * Matrix3_3::Identity();

    Matrix18_18 Q = Matrix18_18::Zero();
    Q.block<3, 3>(0, 0)   = gyro_cov_mat * dt + gyro_bias_cov_mat * (dt3 / 3.0f);
    Q.block<3, 3>(0, 9)   = -gyro_bias_cov_mat * (dt2 / 2.0f);
    Q.block<3, 3>(3, 3)   = accel_cov_mat * (dt3 / 3.0f) + accel_bias_cov_mat * (dt5 / 20.0f);
    Q.block<3, 3>(3, 6)   = accel_cov_mat * (dt2 / 2.0f) + accel_bias_cov_mat * (dt4 / 8.0f);
    Q.block<3, 3>(3, 12)  = -accel_bias_cov_mat * (dt3 / 6.0f);
    Q.block<3, 3>(6, 3)   = Q.block<3, 3>(3, 6);
    Q.block<3, 3>(6, 6)   = accel_cov_mat * dt + accel_bias_cov_mat * (dt3 / 3.0f);
    Q.block<3, 3>(6, 12)  = -accel_bias_cov_mat * (dt2 / 2.0f);
    Q.block<3, 3>(9, 0)   = Q.block<3, 3>(0, 9);
    Q.block<3, 3>(9, 9)   = gyro_bias_cov_mat * dt;
    Q.block<3, 3>(12, 3)  = Q.block<3, 3>(3, 12);
    Q.block<3, 3>(12, 6)  = Q.block<3, 3>(6, 12);
    Q.block<3, 3>(12, 12) = accel_bias_cov_mat * dt;
    Q.block<3, 3>(15, 15) = mag_bias_cov_mat * dt;
    return Q;
}
