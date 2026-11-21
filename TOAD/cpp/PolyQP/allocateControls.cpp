// Transpiled from TOAD/Guidance/PolynomialQP/allocateControls.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

PolyQPStatus allocateControls(const WrenchT &Wrench, const LimitsT &Limits, Mat4Grid &U) {
    (void)Limits; // ASTRA-ONLY: Limits was read only for the vehicle selection
    // Algebraic Control Allocation Matching 6-DoF Plant Dynamics
    const Mat3Grid &F_body = Wrench.F_body;
    const Mat3Grid &Tau_body = Wrench.Tau_body;
    const RowGrid &CGz_vec = Wrench.CGz;
    const int nPts = static_cast<int>(F_body.cols());
    U.setZero(4, nPts);

    for (int k = 1; k <= nPts; ++k) {
        const Vec3 Fk = Col1(F_body, k);
        const Vec3 tauk = Col1(Tau_body, k);
        const Real l_k = mmax(Real(1e-3), V1(CGz_vec, k));

        // Axial force from translational dynamics
        const Real Fz = mmax(Real(1e-3), V1(Fk, 3));

        // ASTRA-ONLY: the vehicle 0 branch below is the only one; the vehicle 1 (TOAD) branch was removed.
        // Vehicle 0 (ASTRA): M_B = [0; 0; -l] x T_B + roll * t_dir
        // Coupling: tau_x = l*Fy + roll*(Fx/T), tau_y = -l*Fx + roll*(Fy/T), tau_z = roll*(Fz/T)
        // Initialize from pure moment arm, then 2-step exact fixed-point refinement for roll * t_dir
        Real Fx = -V1(tauk, 2) / l_k;
        Real Fy = V1(tauk, 1) / l_k;
        Real T_mag = std::sqrt(std::pow(Fx, Real(2)) + std::pow(Fy, Real(2)) + std::pow(Fz, Real(2)));
        Real roll_cmd = V1(tauk, 3) * (T_mag / Fz);
        for (int it = 1; it <= 2; ++it) {
            Fx = (-V1(tauk, 2) + roll_cmd * (Fy / T_mag)) / l_k;
            Fy = (V1(tauk, 1) - roll_cmd * (Fx / T_mag)) / l_k;
            T_mag = std::sqrt(std::pow(Fx, Real(2)) + std::pow(Fy, Real(2)) + std::pow(Fz, Real(2)));
            roll_cmd = V1(tauk, 3) * (T_mag / Fz);
        }

        // Exact trigonometric gimbal inversion for T_B = T * [cos(th)*sin(ph); -sin(th); cos(th)*cos(ph)]
        const Real sin_th = mmax(Real(-1), mmin(Real(1), -Fy / T_mag));
        const Real theta = std::asin(sin_th);
        const Real phi = std::atan2(Fx, Fz);

        // TRANSPILE-FLAG F30: U(4) is a roll torque (N m), consistent with the plant and the LQR model; the controller labels it rad/s^2 and clamps +-4, and evaluateConstraints compares it to RollLimit = 0.9 * MaxRollRate (BLOCKER B1 in ABORT_FIXES.md). Transpiled as-is per MATLAB.
        Col1(U, k) = Vec4(theta, phi, T_mag, roll_cmd);
    }
    return PolyQPStatus::Ok;
}

}
