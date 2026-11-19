// Transpiled from TOAD/Guidance/PolynomialQP/PolyEval.m @ 5ac6bd23
#include "PolyQP.h"
#include "matlab_compat.h"

namespace PolyQP {

PolyQPStatus PolyEval(const QPTrajT &QPTraj, const RowGrid &t, Mat3Grid &Out, int k) {
    // Stateless Evaluation of k-th Derivative of Physical-Time Polynomial
    const RowGrid &t_row = t;
    const int nPts = static_cast<int>(t_row.size());
    const int M = QPTraj.Count;
    const int D = QPTraj.Coef.size(3);
    const int nc = 2 * QPTraj.R;

    Out.setZero(D, nPts);
    if (k >= nc || M == 0) {
        return PolyQPStatus::Ok;
    }

    // Locate active segment index for each query time (favour flight/ramp segments at knot boundaries)
    const RowKnot kt = QPTraj.KeyTimes(R1(1, M + 1));
    const RowSeg T_seg = QPTraj.T(R1(1, M));
    // TRANSPILE-FLAG F13: isfield(QPTraj, 'IsCruise') is always true for a QPTraj from solveAllLegs; transpiled as constant true. Transpiled as-is per MATLAB.
    const bool hasCruise = true;
    static IntGrid idx;
    idx.setConstant(nPts, 1);
    for (int q = 1; q <= nPts; ++q) {
        const Real tq = V1(t_row, q);
        int m_found = 1;
        for (int m = 1; m <= M; ++m) {
            if (tq > V1(kt, m)) {
                m_found = m;
            } else if (tq == V1(kt, m)) {
                const bool isSpecial = V1(QPTraj.IsHold, m) || (hasCruise && V1(QPTraj.IsCruise, m));
                if (!isSpecial) {
                    m_found = m;
                }
            }
        }
        V1(idx, q) = m_found;
    }

    // Local segment time t_loc in [0, T_m]
    static RowGrid t_loc;
    t_loc.resize(nPts);
    for (int q = 1; q <= nPts; ++q) V1(t_loc, q) = V1(t_row, q) - V1(kt, V1(idx, q));
    for (int q = 1; q <= nPts; ++q) {
        if (V1(t_loc, q) < 0) {
            V1(t_loc, q) = 0;
        } else if (V1(t_loc, q) > V1(T_seg, V1(idx, q))) {
            V1(t_loc, q) = V1(T_seg, V1(idx, q));
        }
    }

    // Evaluate k-th derivative in physical time
    static RowGrid tp;
    static RowGrid c_row;
    for (int i = k; i <= nc - 1; ++i) {
        const Real f = factorial(i) / factorial(i - k);
        tp.resize(nPts);
        for (int q = 1; q <= nPts; ++q) V1(tp, q) = std::pow(V1(t_loc, q), Real(i - k));
        for (int d = 1; d <= D; ++d) {
            c_row.resize(nPts);
            for (int q = 1; q <= nPts; ++q) V1(c_row, q) = QPTraj.Coef(i + 1, V1(idx, q), d);
            Row1(Out, d) = Row1(Out, d) + f * (c_row.array() * tp.array()).matrix();
        }
    }
    return PolyQPStatus::Ok;
}

}
