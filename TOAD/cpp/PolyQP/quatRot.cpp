// Transpiled from TOAD/Helper/quatRot.m @ 5ac6bd23
#include "PolyQP.h"
#include "matlab_compat.h"

namespace PolyQP {

Mat3 quatRot(const Vec4 &q) {
    // TRANSPILE-FLAG F36: duplicates the float quatRot in cpp/matlab_helpers.cpp (Real here); the two must stay in sync. Added for FIX-F03.
    const Real q0 = V1(q, 1);
    const Real q1 = V1(q, 2);
    const Real q2 = V1(q, 3);
    const Real q3 = V1(q, 4);
    Mat3 CIB;
    CIB << 1 - 2 * (std::pow(q2, Real(2)) + std::pow(q3, Real(2))), 2 * (q1 * q2 + q0 * q3), 2 * (q1 * q3 - q0 * q2),
        2 * (q1 * q2 - q0 * q3), 1 - 2 * (std::pow(q1, Real(2)) + std::pow(q3, Real(2))), 2 * (q2 * q3 + q0 * q1),
        2 * (q1 * q3 + q0 * q2), 2 * (q2 * q3 - q0 * q1), 1 - 2 * (std::pow(q1, Real(2)) + std::pow(q2, Real(2)));
    return CIB;
}

}
