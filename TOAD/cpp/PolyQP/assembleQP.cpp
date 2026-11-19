// Transpiled from TOAD/Guidance/PolynomialQP/assembleQP.m @ 5ac6bd23
#include "PolyQP.h"
#include "matlab_compat.h"

namespace PolyQP {

PolyQPStatus assembleQP(int r, MatNC &A0, MatNC &Q0, MatNC &H0) {
    // Assemble Normalised-Time QP Matrices (Paper Eq. ??)
    // TRANSPILE-FLAG F14: "Paper Eq. ??" equation number missing in the source comment; preserved verbatim. Transpiled as-is per MATLAB.
    const int nc = 2 * r;
    POLYQP_CAPACITY(r >= 0 && r <= POLYQP_MAX_R && nc <= POLYQP_MAX_NC);
    A0.setZero(nc, nc);
    Q0.setZero(nc, nc);

    // Endpoint derivative map A0 at tau = 0 and tau = 1
    for (int k = 0; k <= r - 1; ++k) {
        M1(A0, k + 1, k + 1) = factorial(k);
        for (int i = k; i <= nc - 1; ++i) {
            M1(A0, r + k + 1, i + 1) = factorial(i) / factorial(i - k);
        }
    }

    // Cost Hessian Q0 for integral of squared r-th derivative over [0, 1]
    for (int i = r; i <= nc - 1; ++i) {
        const Real fi = factorial(i) / factorial(i - r);
        for (int j = r; j <= nc - 1; ++j) {
            const Real fj = factorial(j) / factorial(j - r);
            M1(Q0, i + 1, j + 1) = (fi * fj) / Real(i + j - 2 * r + 1);
        }
    }

    // Exact block inversion of A0 = [D0, 0; A10, A11] to preserve PSD Hessian H0
    // Only the bottom r rows of A0^-1 (corresponding to c_r .. c_{2r-1}) enter Q0
    const MatR A10 = A0(R1(r + 1, nc), R1(1, r));
    const MatR A11 = A0(R1(r + 1, nc), R1(r + 1, nc));
    const MatR D0_inv = A0(R1(1, r), R1(1, r)).diagonal().cwiseInverse().asDiagonal();

    const MatR B1 = Eigen::PartialPivLU<MatR>(A11).solve(MatR::Identity(r, r));
    const MatR B0 = -B1 * (A10 * D0_inv);
    MatRNC B_bot(r, nc);
    B_bot << B0, B1;

    const MatR Q_rr = Q0(R1(r + 1, nc), R1(r + 1, nc));
    H0 = B_bot.transpose() * Q_rr * B_bot;
    const MatNC H0_T = H0.transpose();
    H0 = Real(0.5) * (H0 + H0_T);
    return PolyQPStatus::Ok;
}

}
