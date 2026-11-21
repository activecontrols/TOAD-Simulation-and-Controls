// Transpiled from TOAD/Guidance/PolynomialQP/allocateTimes.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

static void softmaxTimes(const VecSeg &a_log, Real Ttot, Real MinSegTime, RowSeg &T);

PolyQPStatus allocateTimes(const Fix3 &FixLeg, const RowSeg &SegTimesInit, const MatNC &A0, const MatNC &H0, Real MinSegTime, RowSeg &SegTimes, Real &Cost) {
    // Analytic-Gradient Time Redistribution (Paper Sec. V)
    static Coef3 CoefIgnored;
    static Fix3 DerivIgnored;
    SegTimes = SegTimesInit;
    const int M = static_cast<int>(SegTimes.size());
    Real J0;
    POLYQP_TRY(solveLeg(FixLeg, SegTimes, A0, H0, CoefIgnored, DerivIgnored, J0));
    Cost = J0;
    if (M < 2) {
        return PolyQPStatus::Ok;
    }

    const int r = FixLeg.size(2);
    const int nc = 2 * r;
    const Real Ttot = sumv(SegTimes);
    VecSeg a_log = mlog(SegTimes.transpose());
    const int MaxIter = 25;

    // Derivative order vector k_i for segment endpoint vector [0..r-1, 0..r-1]
    // Physical integral Hessian: H_m(T_m) = T_m^(k_i + k_j + 1 - 2r) * H0
    VecNC k_vec(nc);
    for (int q = 1; q <= r; ++q) {
        V1(k_vec, q) = Real(q - 1);
        V1(k_vec, r + q) = Real(q - 1);
    }
    MatNC K_sum(nc, nc);
    for (int j = 1; j <= nc; ++j) {
        for (int i = 1; i <= nc; ++i) M1(K_sum, i, j) = V1(k_vec, i) + V1(k_vec, j);
    }
    const MatNC ExpMat = (Real(1 - 2 * r) + K_sum.array()).matrix();

    for (int it = 1; it <= MaxIter; ++it) {
        // Solve QP at current durations to get optimal endpoint derivatives d*
        RowSeg T_cur;
        softmaxTimes(a_log, Ttot, MinSegTime, T_cur);
        static Fix3 Deriv;
        POLYQP_TRY(solveLeg(FixLeg, T_cur, A0, H0, CoefIgnored, Deriv, J0));

        // Analytic gradient dJ*/dT_m = d*^T * (dH_m/dT_m) * d* (envelope theorem)
        VecSeg grad_T = VecSeg::Zero(M);
        for (int m = 1; m <= M; ++m) {
            const Real Tm = V1(T_cur, m);
            const MatNC dHm_dTm =
                ((ExpMat.array() * (ExpMat.array() - Real(1)).unaryExpr([Tm](Real e) { return std::pow(Tm, e); })) * H0.array()).matrix();
            Real gm = 0;
            for (int d = 1; d <= FixLeg.size(1); ++d) {
                VecNC d_seg;
                reshapeRowSlice(Deriv, d, m, m + 1, d_seg);
                gm = gm + (d_seg.transpose() * dHm_dTm * d_seg).value();
            }
            V1(grad_T, m) = gm;
        }

        // Chain rule to log-durations: dJ/da_m = T_m * (dJ/dT_m - sum_j (T_j/Ttot) * dJ/dT_j)
        const Real mean_grad = sumv(((T_cur.transpose() / Ttot).array() * grad_T.array()).matrix());
        VecSeg g = (T_cur.transpose().array() * (grad_T.array() - mean_grad)).matrix();
        g = (g.array() - mean(g)).matrix();

        if (maxv(mabs(g)) < Real(1e-8) * mmax(J0, Real(1e-12))) {
            break;
        }

        const VecSeg gn = g / maxv(mabs(g));
        bool improved = false;
        Eigen::Matrix<Real, 1, 5> steps;
        steps << Real(0.4), Real(0.2), Real(0.1), Real(0.05), Real(0.025);
        for (int s_idx = 1; s_idx <= static_cast<int>(steps.size()); ++s_idx) {
            const VecSeg an = a_log - V1(steps, s_idx) * gn;
            RowSeg Tn;
            softmaxTimes(an, Ttot, MinSegTime, Tn);
            Real Jn;
            POLYQP_TRY(solveLeg(FixLeg, Tn, A0, H0, CoefIgnored, DerivIgnored, Jn));
            if (Jn < J0 * (1 - Real(1e-6))) {
                a_log = an;
                J0 = Jn;
                improved = true;
                break;
            }
        }

        if (!improved) {
            break;
        }
    }

    softmaxTimes(a_log, Ttot, MinSegTime, SegTimes);
    Cost = J0;
    return PolyQPStatus::Ok;
}

static void softmaxTimes(const VecSeg &a_log, Real Ttot, Real MinSegTime, RowSeg &T) {
    // Map log-durations to positive segment times summing to Ttot
    const VecSeg w = mexp((a_log.array() - maxv(a_log)).matrix());
    T = Ttot * (w.transpose() / sumv(w));
    T = mmax(T, MinSegTime);
    T = Ttot * (T / sumv(T));
}

}
