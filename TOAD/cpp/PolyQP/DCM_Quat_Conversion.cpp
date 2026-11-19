// Transpiled from TOAD/Helper/DCM_Quat_Conversion.m @ 5ac6bd23
// TRANSPILE-FLAG F15: duplicates the float DCM_Quat_Conversion in cpp/matlab_helpers.cpp; the two must stay in sync. Transpiled as-is per MATLAB.
#include "PolyQP.h"
#include "matlab_compat.h"

namespace PolyQP {

// Convert a DCM to a quaternion robustly
Vec4 DCM_Quat_Conversion(const Mat3 &R) {
    // Extract trace
    const Real tr = M1(R, 1, 1) + M1(R, 2, 2) + M1(R, 3, 3);

    // Pre-allocate output
    Vec4 q = Vec4::Zero();

    // CASE 1: Standard Case (Trace > 0)
    // This works for most orientations (small angles, < 90 deg tilt)
    if (tr > 0) {
        const Real S = std::sqrt(tr + Real(1.0)) * 2; // S = 4 * qw
        V1(q, 1) = Real(0.25) * S;
        V1(q, 2) = (M1(R, 3, 2) - M1(R, 2, 3)) / S;
        V1(q, 3) = (M1(R, 1, 3) - M1(R, 3, 1)) / S;
        V1(q, 4) = (M1(R, 2, 1) - M1(R, 1, 2)) / S;

    }
    // CASE 2: Singularity Avoidance (Trace <= 0)
    // We must find the largest diagonal element to avoid dividing by zero.
    else {
        if ((M1(R, 1, 1) > M1(R, 2, 2)) && (M1(R, 1, 1) > M1(R, 3, 3))) {
            // Column 1 (X) is dominant
            const Real S = std::sqrt(Real(1.0) + M1(R, 1, 1) - M1(R, 2, 2) - M1(R, 3, 3)) * 2;
            V1(q, 1) = (M1(R, 3, 2) - M1(R, 2, 3)) / S;
            V1(q, 2) = Real(0.25) * S;
            V1(q, 3) = (M1(R, 1, 2) + M1(R, 2, 1)) / S;
            V1(q, 4) = (M1(R, 1, 3) + M1(R, 3, 1)) / S;

        } else if (M1(R, 2, 2) > M1(R, 3, 3)) {
            // Column 2 (Y) is dominant
            const Real S = std::sqrt(Real(1.0) + M1(R, 2, 2) - M1(R, 1, 1) - M1(R, 3, 3)) * 2;
            V1(q, 1) = (M1(R, 1, 3) - M1(R, 3, 1)) / S;
            V1(q, 2) = (M1(R, 1, 2) + M1(R, 2, 1)) / S;
            V1(q, 3) = Real(0.25) * S;
            V1(q, 4) = (M1(R, 2, 3) + M1(R, 3, 2)) / S;

        } else {
            // Column 3 (Z) is dominant
            const Real S = std::sqrt(Real(1.0) + M1(R, 3, 3) - M1(R, 1, 1) - M1(R, 2, 2)) * 2;
            V1(q, 1) = (M1(R, 2, 1) - M1(R, 1, 2)) / S;
            V1(q, 2) = (M1(R, 1, 3) + M1(R, 3, 1)) / S;
            V1(q, 3) = (M1(R, 2, 3) + M1(R, 3, 2)) / S;
            V1(q, 4) = Real(0.25) * S;
        }
    }

    // Enforce normalization to correct numerical drift
    q = q / mnorm(q);
    return q;
}

}
