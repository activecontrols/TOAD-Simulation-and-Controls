// Transpiled from TOAD/Guidance/PolynomialQP/PolyEval.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

Vector3d PolyEvalPoint(const QPTrajT &QPTraj, double t, int k) {
    // Stateless Evaluation of k-th Derivative of Physical-Time Polynomial
    // STREAM: PolyEval.m's body for one query time; PolyEval (grid) below calls it for each time, so both give the same numbers.
    const int M = QPTraj.Count;
    const int D = QPTraj.Coef.size(3);
    const int nc = 2 * QPTraj.R;
    POLYQP_ASSERT(D == 3);

    Vector3d Out = Vector3d::Zero();
    if (k >= nc || M == 0) {
        return Out;
    }

    // Locate active segment index for each query time (favour flight/ramp segments at knot boundaries)
    // TRANSPILE-FLAG F13: isfield(QPTraj, 'IsCruise') is always true for a QPTraj from solveAllLegs; transpiled as constant true. Transpiled as-is per MATLAB.
    const bool hasCruise = true;
    const double tq = t;
    int m_found = 1;
    for (int m = 1; m <= M; ++m) {
        if (tq > V1(QPTraj.KeyTimes, m)) {
            m_found = m;
        } else if (tq == V1(QPTraj.KeyTimes, m)) {
            const bool isSpecial = V1(QPTraj.IsHold, m) || (hasCruise && V1(QPTraj.IsCruise, m));
            if (!isSpecial) {
                m_found = m;
            }
        }
    }
    const int idx = m_found;

    // Local segment time t_loc in [0, T_m]
    double t_loc = tq - V1(QPTraj.KeyTimes, idx);
    if (t_loc < 0) {
        t_loc = 0;
    } else if (t_loc > V1(QPTraj.T, idx)) {
        t_loc = V1(QPTraj.T, idx);
    }

    // Evaluate k-th derivative in physical time
    for (int i = k; i <= nc - 1; ++i) {
        const double f = factorial(i) / factorial(i - k);
        const double tp = std::pow(t_loc, double(i - k));
        for (int d = 1; d <= D; ++d) {
            const double c_row = QPTraj.Coef(i + 1, idx, d);
            V1(Out, d) = V1(Out, d) + f * (c_row * tp);
        }
    }
    return Out;
}

PolyQPStatus PolyEval(const QPTrajT &QPTraj, const Matrix1_Nd &t, int nPts, Matrix3_Nd &Out, int k) {
    // STREAM: the grid form (kept for verification): PolyEvalPoint at each of the nPts query times.
    POLYQP_CAPACITY(nPts <= POLYQP_MAX_NODES);
    Out.setZero();
    for (int q = 1; q <= nPts; ++q) {
        Col1(Out, q) = PolyEvalPoint(QPTraj, V1(t, q), k);
    }
    return PolyQPStatus::Ok;
}

}
