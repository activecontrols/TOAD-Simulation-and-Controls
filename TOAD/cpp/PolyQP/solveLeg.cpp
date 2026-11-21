// Transpiled from TOAD/Guidance/PolynomialQP/solveLeg.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

// STATIC: the rcond guard is evaluated on the leading nP x nP block of R_PP_eq with a fixed-size LU of exactly that size
// (N = 1..POLYQP_MAX_FREE, chosen at run time). The identity padding leaves the exact 1-norm condition number unchanged, but
// Eigen's rcond is an estimate (Hager / Higham), and on the padded matrix it differed by up to 3x on 27 of 1874 fuzzed axis
// solves (verify_polyqp, Phase 4). Estimating on the unpadded block keeps the guard identical to the variable-size solver.
template <int N> struct RcondLeading {
    static double eval(const MatrixF_Fd &A, int n) {
        if (n == N) {
            const Eigen::Matrix<double, N, N> B = A.template topLeftCorner<N, N>();
            return Eigen::PartialPivLU<Eigen::Matrix<double, N, N>>(B).rcond();
        }
        return RcondLeading<N - 1>::eval(A, n);
    }
};
template <> struct RcondLeading<0> {
    static double eval(const MatrixF_Fd &, int) { return 0.0; }
};

PolyQPStatus solveLeg(const Array3_5_Kd &FixLeg, const Array3_5_Kb &FixLegFixed, const Matrix1_Sd &SegTimes, const Matrix10_10d &A0, const Matrix10_10d &H0,
                      Array10_S_3d &CoefPhys, Array3_5_Kd &Deriv, double &Cost) {
    // Closed-Form Minimum-Derivative QP Solve (Paper Eq. ??)
    // TRANSPILE-FLAG F14: "Paper Eq. ??" equation number missing in the source comment; preserved verbatim. Transpiled as-is per MATLAB.
    const int D = FixLeg.size(1);
    const int r = FixLeg.size(2);
    const int K = FixLeg.size(3);
    const int M = K - 1;
    const int nc = 2 * r;
    const int nv = r * K;
    POLYQP_CAPACITY(r == POLYQP_R && nv <= POLYQP_MAX_NV && M <= POLYQP_MAX_SEG);

    // STATIC: the buffers live in the shared workspace (PolyQP_types.h); the solver uses their leading nv x nv / nP x nP blocks.
    PolyQPWorkspace::SolveLegT &W = workspace().solveLeg;

    // Assemble global endpoint-derivative Hessian Rmat over leg knots
    // Physical action integral weight T_m^(1-2r) for squared r-th derivative over [0, T_m]
    MatrixNV_NVd &Rmat = W.Rmat;
    Rmat.setZero();
    Array10_10_Sd &Hseg = W.Hseg;
    Hseg.setConstant(nc, nc, M, double(0)); // FIX-F07: per-segment Hessians for the cancellation-free reduced assembly
    for (int m = 1; m <= M; ++m) {
        const double Tm = V1(SegTimes, m);
        Vector10d s;
        for (int q = 1; q <= r; ++q) {
            V1(s, q) = std::pow(Tm, double(q - 1));
            V1(s, r + q) = std::pow(Tm, double(q - 1));
        }
        const Matrix10_10d s_sT = s * s.transpose();
        const Matrix10_10d Hm = ((std::pow(Tm, double(1 - 2 * r)) * s_sT).array() * H0.array()).matrix();
        const auto ix = R1((m - 1) * r + 1, (m - 1) * r + nc);
        Rmat(ix, ix) = Rmat(ix, ix) + Hm;
        for (int j = 1; j <= nc; ++j) {
            for (int i = 1; i <= nc; ++i) Hseg(i, j, m) = M1(Hm, i, j); // FIX-F07: per-segment Hessians for the cancellation-free reduced assembly
        }
    }

    // Solve for free derivatives dP per axis: dP* = -R_PP \ (R_FP' * dF)
    Deriv.setConstant(D, r, K, double(0));
    Cost = 0;
    for (int d = 1; d <= D; ++d) {
        VectorNVd &fx = W.fx;
        reshapeRowSlice(FixLeg, d, 1, K, fx);
        VectorNVb &fxFixed = W.fxFixed; // NO-NAN: ~isnan(fx)
        reshapeRowSlice(FixLegFixed, d, 1, K, fxFixed);
        const FindResult<POLYQP_MAX_NV> idxF = find<POLYQP_MAX_NV>(fxFixed.head(nv).array());
        const FindResult<POLYQP_MAX_NV> idxP = find<POLYQP_MAX_NV>(!fxFixed.head(nv).array());
        const int nF = idxF.numel();
        const int nP = idxP.numel();

        if (idxF.numel() < r) {
            POLYQP_ERROR(IllPosedQP, "Axis %d has %d fixed BCs; at least r = %d required to pin null space.", d, idxF.numel(), r);
        }
        POLYQP_CAPACITY(nP <= POLYQP_MAX_FREE); // STATIC: R_PP is at most POLYQP_MAX_FREE x POLYQP_MAX_FREE

        VectorNVd &dv = W.dv;
        dv = fx;
        Array1_Sb isShiftSeg = Array1_Sb::Constant(false); // FIX-F07: segments of this axis whose cost is invariant to a common shift of both positions
        MatrixNV_NVd &Rns = W.Rns;
        Rns.setZero(); // FIX-F07: Hessian of the segments that are not shift-invariant
        if (!idxP.isempty()) {
            // STATIC: R_PP, R_PF, Tz and the equilibrated system are fixed POLYQP_MAX_FREE-sized; the leading nP rows / columns are used.
            MatrixF_Fd &R_PP = W.R_PP;
            R_PP.topLeftCorner(nP, nP) = double(0.5) * (Rmat(idxP.zb(), idxP.zb()) + Rmat(idxP.zb(), idxP.zb()).transpose());
            MatrixF_NVd &R_PF = W.R_PF;
            R_PF.topLeftCorner(nP, nF) = Rmat(idxP.zb(), idxF.zb());
            VectorNVd &dF = W.dF;
            dF.head(nF) = fx(idxF.zb());

            // FIX-F07: a segment with derivatives 1..r-1 fixed at both knots and both positions free has a cost invariant to a common shift; reparametrise p_{m+1} = p_m + Delta_m (d_P = Tz * z)
            MatrixF_Fd &Tz = W.Tz;
            Tz.setIdentity();
            for (int m = 1; m <= M; ++m) {
                const int iA = (m - 1) * r + 1;
                const int iB = m * r + 1;
                if (!V1(fxFixed, iA) && !V1(fxFixed, iB) && fxFixed(R1(iA + 1, iA + r - 1)).all() && fxFixed(R1(iB + 1, iB + r - 1)).all()) { // NO-NAN: isnan(fx(i)) is !fxFixed(i)
                    const int jA = findFirst(idxP, iA);
                    const int jB = findFirst(idxP, iB);
                    Tz.row(jB - 1) = Tz.row(jA - 1);
                    M1(Tz, jB, jB) = double(1);
                    V1(isShiftSeg, m) = true;
                }
            }
            if (isShiftSeg.head(M).any()) { // FIX-F07: reduced system Tz' * R_PP * Tz, assembled so that the common-shift direction cancels exactly
                for (int m = 1; m <= M; ++m) {
                    if (!V1(isShiftSeg, m)) {
                        for (int j = 1; j <= nc; ++j) {
                            for (int i = 1; i <= nc; ++i) M1(Rns, (m - 1) * r + i, (m - 1) * r + j) = M1(Rns, (m - 1) * r + i, (m - 1) * r + j) + Hseg(i, j, m);
                        }
                    }
                }
                MatrixF_Fd &Rns_PP = W.Rns_PP;
                Rns_PP.topLeftCorner(nP, nP) = double(0.5) * (Rns(idxP.zb(), idxP.zb()) + Rns(idxP.zb(), idxP.zb()).transpose());
                MatrixF_Fd &TzT_Rns_PP = W.TzT_Rns_PP;
                TzT_Rns_PP.topLeftCorner(nP, nP).noalias() = Tz.topLeftCorner(nP, nP).transpose() * Rns_PP.topLeftCorner(nP, nP);
                R_PP.topLeftCorner(nP, nP).noalias() = TzT_Rns_PP.topLeftCorner(nP, nP) * Tz.topLeftCorner(nP, nP);
                MatrixF_Fd &R_PP_T = W.R_PP_T;
                R_PP_T.topLeftCorner(nP, nP) = R_PP.topLeftCorner(nP, nP).transpose();
                R_PP.topLeftCorner(nP, nP) = double(0.5) * (R_PP.topLeftCorner(nP, nP) + R_PP_T.topLeftCorner(nP, nP));
                MatrixF_NVd &Rns_PF = W.Rns_PF;
                Rns_PF.topLeftCorner(nP, nF) = Rns(idxP.zb(), idxF.zb());
                R_PF.topLeftCorner(nP, nF).noalias() = Tz.topLeftCorner(nP, nP).transpose() * Rns_PF.topLeftCorner(nP, nF);
                // FIX-F07: a shift-invariant segment enters only through Delta_m (its p_{m+1} row and column, local index r + 1)
                for (int m = 1; m <= M; ++m) {
                    if (V1(isShiftSeg, m)) {
                        const int jB = findFirst(idxP, m * r + 1);
                        M1(R_PP, jB, jB) = M1(R_PP, jB, jB) + Hseg(r + 1, r + 1, m);
                        for (int q = 1; q <= idxF.numel(); ++q) {
                            const int iL = idxF(q) - (m - 1) * r;
                            if (iL >= 1 && iL <= nc) {
                                M1(R_PF, jB, q) = M1(R_PF, jB, q) + Hseg(r + 1, iL, m);
                            }
                        }
                    }
                }
            }

            // Jacobi diagonal equilibration: D_inv = 1 ./ sqrt(diag(R_PP))
            VectorFd &d_diag = W.d_diag;
            d_diag.head(nP) = R_PP.diagonal().head(nP);
            // NO-NAN: the MATLAB guard also rejects a non-finite diagonal; with real inputs it is finite, so only the numerical d_diag <= 0 test remains
            if ((d_diag.head(nP).array() <= 0).any()) {
                POLYQP_ERROR(SingularRPP, "Axis %d R_PP has non-positive diagonal (min = %.2e).", d, static_cast<double>(minv(d_diag.head(nP))));
            }
            VectorFd &s_jac = W.s_jac;
            s_jac.head(nP) = msqrt(d_diag.head(nP)).cwiseInverse();
            MatrixF_Fd &s_jac_s_jacT = W.s_jac_s_jacT;
            s_jac_s_jacT.topLeftCorner(nP, nP).noalias() = s_jac.head(nP) * s_jac.head(nP).transpose();
            // STATIC: R_PP_eq is solved as a fixed POLYQP_MAX_FREE x POLYQP_MAX_FREE system: the nP x nP equilibrated matrix padded
            // STATIC: with an identity block and a zero right-hand side. The system is block-diagonal, so the padding cannot change
            // STATIC: the first nP unknowns; and R_PP_eq has a unit diagonal, so its largest eigenvalue is >= 1 and its smallest <= 1:
            // STATIC: the added unit eigenvalues leave its condition number (and the rcond guard) unchanged (verify_polyqp checks both).
            MatrixF_Fd &R_PP_eq = W.R_PP_eq;
            R_PP_eq.setIdentity();
            R_PP_eq.topLeftCorner(nP, nP) = (s_jac_s_jacT.topLeftCorner(nP, nP).array() * R_PP.topLeftCorner(nP, nP).array()).matrix();
            MatrixF_Fd &R_PP_eq_T = W.R_PP_eq_T;
            R_PP_eq_T.topLeftCorner(nP, nP) = R_PP_eq.topLeftCorner(nP, nP).transpose();
            R_PP_eq.topLeftCorner(nP, nP) = double(0.5) * (R_PP_eq.topLeftCorner(nP, nP) + R_PP_eq_T.topLeftCorner(nP, nP));
            VectorFd &R_PF_dF = W.R_PF_dF;
            R_PF_dF.head(nP).noalias() = R_PF.topLeftCorner(nP, nF) * dF.head(nF);
            VectorFd &rhs_eq = W.rhs_eq;
            rhs_eq.setZero();
            rhs_eq.head(nP) = (s_jac.head(nP).array() * R_PF_dF.head(nP).array()).matrix();

            const double rc = RcondLeading<POLYQP_MAX_FREE>::eval(R_PP_eq, nP); // STATIC: rcond of the unpadded nP x nP block
#if defined(POLYQP_ERROR_MESSAGES)
            // STATIC: verification only: the guard's rcond, the padded matrix's own estimate, and the system, per axis
            lastRcond().rc[d - 1] = rc;
            lastRcond().rcPadded[d - 1] = Eigen::PartialPivLU<MatrixF_Fd>(R_PP_eq).rcond();
            lastRcond().nP[d - 1] = nP;
            lastRcond().R_PP_eq[d - 1] = R_PP_eq;
#endif
            // TRANSPILE-FLAG F07 (fixed by FIX-F07): lateral R_PP was near-singular for unbalanced segment times (cond ~ 4e13 at (10, 0.6) s); SingularRPP could trip even in double.
            // NO-NAN: MATLAB tests isnan(rc) || rc < 1e-13. Eigen's rcond estimate of an exactly singular R_PP_eq is NaN (0 / 0 in the
            // NO-NAN: inverse-norm estimate) even for real inputs, so the threshold is written as "rc is not >= 1e-13": the same guard
            // NO-NAN: without an isnan call. A literal rc < 1e-13 let 8 exactly singular fuzzed legs through (verify_polyqp, Phase 3).
            if (!(rc >= double(1e-13))) {
                POLYQP_ERROR(SingularRPP, "Axis %d free-derivative Hessian R_PP is singular (rcond = %.2e).", d, static_cast<double>(rc));
            }
            Eigen::LLT<MatrixF_Fd> &llt = W.llt;
            llt.compute(R_PP_eq);
            VectorFd &R_PP_eq_rhs_eq = W.R_PP_eq_rhs_eq;
            if (llt.info() == Eigen::Success) {
                R_PP_eq_rhs_eq = llt.solve(rhs_eq);
            } else {
                Eigen::PartialPivLU<MatrixF_Fd> &lu = W.lu;
                lu.compute(R_PP_eq);
                R_PP_eq_rhs_eq = lu.solve(rhs_eq);
            }
            dv(idxP.zb()) = ((-R_PP_eq_rhs_eq.head(nP)).array() * s_jac.head(nP).array()).matrix();
            if (isShiftSeg.head(M).any()) { // FIX-F07: map the reduced solution back, d_P = Tz * z
                VectorFd &z = W.z;
                z.head(nP) = dv(idxP.zb());
                VectorFd &Tz_z = W.Tz_z;
                Tz_z.head(nP).noalias() = Tz.topLeftCorner(nP, nP) * z.head(nP);
                dv(idxP.zb()) = Tz_z.head(nP);
            }
        }

        assignRowSlice(Deriv, d, dv);
        if (isShiftSeg.head(M).any()) { // FIX-F07: cost without the cancelling common-shift terms (a shift-invariant segment is evaluated at p_m = 0, p_{m+1} = Delta_m)
            Cost = Cost + (dv.head(nv).transpose() * Rns.topLeftCorner(nv, nv) * dv.head(nv)).value();
            for (int m = 1; m <= M; ++m) {
                if (V1(isShiftSeg, m)) {
                    Vector10d xm = dv(R1((m - 1) * r + 1, (m - 1) * r + nc));
                    V1(xm, r + 1) = V1(xm, r + 1) - V1(xm, 1);
                    V1(xm, 1) = double(0);
                    Matrix10_10d Hm_m;
                    for (int j = 1; j <= nc; ++j) {
                        for (int i = 1; i <= nc; ++i) M1(Hm_m, i, j) = Hseg(i, j, m);
                    }
                    Cost = Cost + (xm.transpose() * Hm_m * xm).value();
                }
            }
        } else {
            Cost = Cost + (dv.head(nv).transpose() * Rmat.topLeftCorner(nv, nv) * dv.head(nv)).value();
        }
    }

    // Recover physical-time polynomial coefficients c_i,phys = c_i,tau / T_m^i
    const Matrix5_5d D0_inv = A0(R1(1, r), R1(1, r)).diagonal().cwiseInverse().asDiagonal();
    const Matrix5_5d A10 = A0(R1(r + 1, nc), R1(1, r));
    const Matrix5_5d A11 = A0(R1(r + 1, nc), R1(r + 1, nc));

    CoefPhys.setConstant(nc, M, D, double(0));
    for (int m = 1; m <= M; ++m) {
        const double Tm = V1(SegTimes, m);
        Vector10d s;
        for (int q = 1; q <= r; ++q) {
            V1(s, q) = std::pow(Tm, double(q - 1));
            V1(s, r + q) = std::pow(Tm, double(q - 1));
        }
        Vector10d scalePhys;
        for (int q = 1; q <= nc; ++q) V1(scalePhys, q) = std::pow(Tm, double(q - 1));
        for (int d = 1; d <= D; ++d) {
            Vector10d d_seg;
            reshapeRowSlice(Deriv, d, m, m + 1, d_seg);
            d_seg = (s.array() * d_seg.array()).matrix();
            const Vector5d c_top = D0_inv * d_seg(R1(1, r));
            const Vector5d c_bot = Eigen::PartialPivLU<Matrix5_5d>(A11).solve(d_seg(R1(r + 1, nc)) - A10 * c_top);
            Vector10d c_tau;
            c_tau << c_top, c_bot;
            for (int q = 1; q <= nc; ++q) CoefPhys(q, m, d) = V1(c_tau, q) / V1(scalePhys, q);
        }
    }
    return PolyQPStatus::Ok;
}

}
