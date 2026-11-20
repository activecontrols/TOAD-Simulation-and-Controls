// Algebraic Control Allocation for 6DoF. Exact in P-space
#include "AbortCode.h"
#include <cmath>

Vector4 allocateControlsPoint(const AbortWrenchPoint &Wrench, float rTB) {
    Vector3 Fk = Wrench.F_body;
    Vector3 tauk = Wrench.Tau_body;
    float l_k = fmaxf(1e-3f, rTB);

    // Axial force from translational dynamics
    float Fz = fmaxf(1e-3f, Fk(2));

    // Roll allocation
    float Fx = -tauk(1) / l_k;
    float Fy = tauk(0) / l_k;
    float T_mag = sqrtf(Fx * Fx + Fy * Fy + Fz * Fz);
    float roll_cmd = tauk(2) * (T_mag / Fz);
    for (int it = 0; it < 2; ++it) {
        Fx = (-tauk(1) + roll_cmd * (Fy / T_mag)) / l_k;
        Fy = (tauk(0) - roll_cmd * (Fx / T_mag)) / l_k;
        T_mag = sqrtf(Fx * Fx + Fy * Fy + Fz * Fz);
        roll_cmd = tauk(2) * (T_mag / Fz);
    }

    // Gimbal inversion
    float sin_th = fmaxf(-1.0f, fminf(1.0f, -Fy / T_mag));
    float theta = asinf(sin_th);
    float phi = atan2f(Fx, Fz);

    return Vector4(theta, phi, T_mag, roll_cmd);
}
