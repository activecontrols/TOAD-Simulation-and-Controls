// Transpiled from TOAD/Guidance/PolynomialQP/allocateControls.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

Vector4d allocateControlsPoint(const Vector3d &Fk, const Vector3d &tauk, double CGz, const LimitsT &Limits) {
    (void)Limits; // ASTRA-ONLY: Limits was read only for the vehicle selection
    // Algebraic Control Allocation Matching 6-DoF Plant Dynamics
    // STREAM: allocateControls.m's loop body for one sample; allocateControls (grid) below calls it for each sample.
    const double l_k = std::max(double(1e-3), CGz);

    // Axial force from translational dynamics
    const double Fz = std::max(double(1e-3), V1(Fk, 3));

    // ASTRA-ONLY: the vehicle 0 branch below is the only one; the vehicle 1 (TOAD) branch was removed.
    // Vehicle 0 (ASTRA): M_B = [0; 0; -l] x T_B + roll * t_dir
    // Coupling: tau_x = l*Fy + roll*(Fx/T), tau_y = -l*Fx + roll*(Fy/T), tau_z = roll*(Fz/T)
    // Initialize from pure moment arm, then 2-step exact fixed-point refinement for roll * t_dir
    double Fx = -V1(tauk, 2) / l_k;
    double Fy = V1(tauk, 1) / l_k;
    double T_mag = std::sqrt(std::pow(Fx, double(2)) + std::pow(Fy, double(2)) + std::pow(Fz, double(2)));
    double roll_cmd = V1(tauk, 3) * (T_mag / Fz);
    for (int it = 1; it <= 2; ++it) {
        Fx = (-V1(tauk, 2) + roll_cmd * (Fy / T_mag)) / l_k;
        Fy = (V1(tauk, 1) - roll_cmd * (Fx / T_mag)) / l_k;
        T_mag = std::sqrt(std::pow(Fx, double(2)) + std::pow(Fy, double(2)) + std::pow(Fz, double(2)));
        roll_cmd = V1(tauk, 3) * (T_mag / Fz);
    }

    // Exact trigonometric gimbal inversion for T_B = T * [cos(th)*sin(ph); -sin(th); cos(th)*cos(ph)]
    const double sin_th = std::max(double(-1), std::min(double(1), -Fy / T_mag));
    const double theta = std::asin(sin_th);
    const double phi = std::atan2(Fx, Fz);

    // TRANSPILE-FLAG F30: U(4) is a roll torque (N m), consistent with the plant and the LQR model; the controller labels it rad/s^2 and clamps +-4, and evaluateConstraints compares it to RollLimit = 0.9 * MaxRollRate (BLOCKER B1 in ABORT_FIXES.md). Transpiled as-is per MATLAB.
    return Vector4d(theta, phi, T_mag, roll_cmd);
}

PolyQPStatus allocateControls(const WrenchT &Wrench, int nPts, const LimitsT &Limits, Matrix4_Nd &U) {
    // STREAM: the grid form (kept for verification): allocateControlsPoint at each of the nPts samples.
    POLYQP_CAPACITY(nPts <= POLYQP_MAX_NODES);
    U.setZero();
    for (int k = 1; k <= nPts; ++k) {
        Col1(U, k) = allocateControlsPoint(Col1(Wrench.F_body, k), Col1(Wrench.Tau_body, k), V1(Wrench.CGz, k), Limits);
    }
    return PolyQPStatus::Ok;
}

}
