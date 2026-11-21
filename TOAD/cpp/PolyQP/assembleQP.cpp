// Transpiled from TOAD/Guidance/PolynomialQP/assembleQP.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

PolyQPStatus assembleQP(Matrix10_10d &A0, Matrix10_10d &Q0, Matrix10_10d &H0) {
    // Assemble Normalised-Time QP Matrices (Paper Eq. ??)
    // TRANSPILE-FLAG F14: "Paper Eq. ??" equation number missing in the source comment; preserved verbatim. Transpiled as-is per MATLAB.
    const int r = POLYQP_R; // ASTRA-ONLY: the polynomial order is hardcoded (r = 5)
    const int nc = 2 * r;
    A0.setZero();
    Q0.setZero();

    // Endpoint derivative map A0 at tau = 0 and tau = 1
    for (int k = 0; k <= r - 1; ++k) {
        M1(A0, k + 1, k + 1) = factorial(k);
        for (int i = k; i <= nc - 1; ++i) {
            M1(A0, r + k + 1, i + 1) = factorial(i) / factorial(i - k);
        }
    }

    // Cost Hessian Q0 for integral of squared r-th derivative over [0, 1]
    for (int i = r; i <= nc - 1; ++i) {
        const double fi = factorial(i) / factorial(i - r);
        for (int j = r; j <= nc - 1; ++j) {
            const double fj = factorial(j) / factorial(j - r);
            M1(Q0, i + 1, j + 1) = (fi * fj) / double(i + j - 2 * r + 1);
        }
    }

    // Exact block inversion of A0 = [D0, 0; A10, A11] to preserve PSD Hessian H0
    // Only the bottom r rows of A0^-1 (corresponding to c_r .. c_{2r-1}) enter Q0
    const Matrix5_5d A10 = A0(R1(r + 1, nc), R1(1, r));
    const Matrix5_5d A11 = A0(R1(r + 1, nc), R1(r + 1, nc));
    const Matrix5_5d D0_inv = A0(R1(1, r), R1(1, r)).diagonal().cwiseInverse().asDiagonal();

    const Matrix5_5d B1 = Eigen::PartialPivLU<Matrix5_5d>(A11).solve(Matrix5_5d::Identity());
    const Matrix5_5d B0 = -B1 * (A10 * D0_inv);
    Matrix5_10d B_bot;
    B_bot << B0, B1;

    const Matrix5_5d Q_rr = Q0(R1(r + 1, nc), R1(r + 1, nc));
    H0 = B_bot.transpose() * Q_rr * B_bot;
    const Matrix10_10d H0_T = H0.transpose();
    H0 = double(0.5) * (H0 + H0_T);
    return PolyQPStatus::Ok;
}

}
