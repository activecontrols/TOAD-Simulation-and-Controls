// Transpiled from TOAD/Guidance/PolynomialQP/allocateTimes.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

static void softmaxTimes(const VectorSd &a_log, int M, double Ttot, double MinSegTime, Matrix1_Sd &T);

PolyQPStatus allocateTimes(const Array3_5_Kd &FixLeg, const Array3_5_Kb &FixLegFixed, const Matrix1_Sd &SegTimesInit, const Matrix10_10d &A0, const Matrix10_10d &H0,
                           double MinSegTime, Matrix1_Sd &SegTimes, double &Cost) {
    // Analytic-Gradient Time Redistribution (Paper Sec. V)
    PolyQPWorkspace::AllocateTimesT &W = workspace().allocateTimes;
    Array10_S_3d &CoefIgnored = W.CoefIgnored;
    Array3_5_Kd &DerivIgnored = W.DerivIgnored;
    SegTimes = SegTimesInit;
    const int M = FixLeg.size(3) - 1; // STATIC: numel(SegTimes); the fixed-size row holds the leg's M segment times first
    double J0;
    POLYQP_TRY(solveLeg(FixLeg, FixLegFixed, SegTimes, A0, H0, CoefIgnored, DerivIgnored, J0));
    Cost = J0;
    if (M < 2) {
        return PolyQPStatus::Ok;
    }

    const int r = FixLeg.size(2);
    const int nc = 2 * r;
    const double Ttot = sumv(SegTimes.head(M));
    VectorSd a_log = VectorSd::Zero();
    a_log.head(M) = mlog(SegTimes.head(M).transpose());
    const int MaxIter = 25;

    // Derivative order vector k_i for segment endpoint vector [0..r-1, 0..r-1]
    // Physical integral Hessian: H_m(T_m) = T_m^(k_i + k_j + 1 - 2r) * H0
    Vector10d k_vec;
    for (int q = 1; q <= r; ++q) {
        V1(k_vec, q) = double(q - 1);
        V1(k_vec, r + q) = double(q - 1);
    }
    Matrix10_10d K_sum;
    for (int j = 1; j <= nc; ++j) {
        for (int i = 1; i <= nc; ++i) M1(K_sum, i, j) = V1(k_vec, i) + V1(k_vec, j);
    }
    const Matrix10_10d ExpMat = (double(1 - 2 * r) + K_sum.array()).matrix();

    for (int it = 1; it <= MaxIter; ++it) {
        // Solve QP at current durations to get optimal endpoint derivatives d*
        Matrix1_Sd T_cur = Matrix1_Sd::Zero();
        softmaxTimes(a_log, M, Ttot, MinSegTime, T_cur);
        Array3_5_Kd &Deriv = W.Deriv;
        POLYQP_TRY(solveLeg(FixLeg, FixLegFixed, T_cur, A0, H0, CoefIgnored, Deriv, J0));

        // Analytic gradient dJ*/dT_m = d*^T * (dH_m/dT_m) * d* (envelope theorem)
        VectorSd grad_T = VectorSd::Zero();
        for (int m = 1; m <= M; ++m) {
            const double Tm = V1(T_cur, m);
            const Matrix10_10d dHm_dTm =
                ((ExpMat.array() * (ExpMat.array() - double(1)).unaryExpr([Tm](double e) { return std::pow(Tm, e); })) * H0.array()).matrix();
            double gm = 0;
            for (int d = 1; d <= FixLeg.size(1); ++d) {
                Vector10d d_seg;
                reshapeRowSlice(Deriv, d, m, m + 1, d_seg);
                gm = gm + (d_seg.transpose() * dHm_dTm * d_seg).value();
            }
            V1(grad_T, m) = gm;
        }

        // Chain rule to log-durations: dJ/da_m = T_m * (dJ/dT_m - sum_j (T_j/Ttot) * dJ/dT_j)
        const double mean_grad = sumv(((T_cur.head(M).transpose() / Ttot).array() * grad_T.head(M).array()).matrix());
        VectorSd g = VectorSd::Zero();
        g.head(M) = (T_cur.head(M).transpose().array() * (grad_T.head(M).array() - mean_grad)).matrix();
        g.head(M) = (g.head(M).array() - mean(g.head(M))).matrix();

        if (maxv(mabs(g.head(M))) < double(1e-8) * std::max(J0, double(1e-12))) {
            break;
        }

        const VectorSd gn = g / maxv(mabs(g.head(M)));
        bool improved = false;
        Eigen::Matrix<double, 1, 5> steps;
        steps << double(0.4), double(0.2), double(0.1), double(0.05), double(0.025);
        for (int s_idx = 1; s_idx <= static_cast<int>(steps.size()); ++s_idx) {
            const VectorSd an = a_log - V1(steps, s_idx) * gn;
            Matrix1_Sd Tn = Matrix1_Sd::Zero();
            softmaxTimes(an, M, Ttot, MinSegTime, Tn);
            double Jn;
            POLYQP_TRY(solveLeg(FixLeg, FixLegFixed, Tn, A0, H0, CoefIgnored, DerivIgnored, Jn));
            if (Jn < J0 * (1 - double(1e-6))) {
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

    softmaxTimes(a_log, M, Ttot, MinSegTime, SegTimes);
    Cost = J0;
    return PolyQPStatus::Ok;
}

static void softmaxTimes(const VectorSd &a_log, int M, double Ttot, double MinSegTime, Matrix1_Sd &T) {
    // Map log-durations to positive segment times summing to Ttot
    // STATIC: only the leading M entries are segment times
    VectorSd w = VectorSd::Zero();
    w.head(M) = mexp((a_log.head(M).array() - maxv(a_log.head(M))).matrix());
    T.head(M) = Ttot * (w.head(M).transpose() / sumv(w.head(M)));
    T.head(M) = T.head(M).cwiseMax(MinSegTime);
    T.head(M) = Ttot * (T.head(M) / sumv(T.head(M)));
}

}
