#include "matlab_funcs.h"

Matrix18_18 Q_gen() {
    float gyro_cov = 1e-5;
    float gyro_bias_cov = 1e-12;
    float accel_proc_cov = 5e-3;
    float accel_bias_cov = 1e-14;
    float mag_proc_cov = 1e-1;
    float mag_bias_cov = 1e-15; // * 3000;

    // Init constants
    float dt = 0.001;

    // body vectors, defined from flight computer notation
    Vector3 f_global_body = (Vector3() << 0, 0, -9.81).finished(); // Expected gravity vector from 0 (standing upright)
    Vector3 m_global_body = (Vector3() << 1, 0, 0).finished();     // Expected normalized field vector from 0 (standing upright)

    Matrix3_3 gyro_cov_mat = gyro_cov * Matrix3_3::Identity();
    Matrix3_3 gyro_bias_cov_mat = gyro_bias_cov * Matrix3_3::Identity();
    Matrix3_3 accel_cov_mat = accel_proc_cov * Matrix3_3::Identity();
    Matrix3_3 accel_bias_cov_mat = accel_bias_cov * Matrix3_3::Identity();
    Matrix3_3 mag_cov_mat = mag_proc_cov * Matrix3_3::Identity();
    Matrix3_3 mag_bias_cov_mat = mag_bias_cov * Matrix3_3::Identity();

    Matrix18_18 Q = Matrix18_18::Zero();
    Q.block<3, 3>(0, 0) = gyro_cov_mat * dt + gyro_bias_cov_mat * (pow(dt, 3)) / 3.0;
    Q.block<3, 3>(0, 9) = -gyro_bias_cov_mat * (pow(dt, 2)) / 2.0;
    Q.block<3, 3>(3, 3) = accel_cov_mat * (pow(dt, 3)) / 3.0 + accel_bias_cov_mat * (pow(dt, 5)) / 20.0;
    Q.block<3, 3>(3, 6) = accel_cov_mat * (pow(dt, 2)) / 2.0 + accel_bias_cov_mat * (pow(dt, 4)) / 8.0;
    Q.block<3, 3>(3, 12) = -accel_bias_cov_mat * (pow(dt, 3)) / 6.0;
    Q.block<3, 3>(6, 3) = accel_cov_mat * (pow(dt, 2)) / 2.0 + accel_bias_cov_mat * (pow(dt, 4)) / 8.0;
    Q.block<3, 3>(6, 6) = accel_cov_mat * dt + accel_bias_cov_mat * (pow(dt, 3)) / 3.0;
    Q.block<3, 3>(6, 12) = -accel_bias_cov_mat * (pow(dt, 2)) / 2.0;
    Q.block<3, 3>(9, 0) = -gyro_bias_cov_mat * (pow(dt, 2)) / 2.0;
    Q.block<3, 3>(9, 9) = gyro_bias_cov_mat * dt;
    Q.block<3, 3>(12, 3) = -accel_bias_cov_mat * (pow(dt, 2)) / 2.0;
    Q.block<3, 3>(12, 6) = -accel_bias_cov_mat * (pow(dt, 3)) / 6.0;
    Q.block<3, 3>(12, 12) = accel_bias_cov_mat * dt;
    Q.block<3, 3>(15, 15) = mag_bias_cov_mat * dt;
    return Q;
}