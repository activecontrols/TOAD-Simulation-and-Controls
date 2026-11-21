// Transpiled from TOAD/Guidance/PolynomialQP/solveLeg.m @ 5ac6bd23
#include "PolyQP.h"

namespace PolyQP {

PolyQPStatus solveLeg(const Fix3 &FixLeg, const Bool3 &FixLegFixed, const RowSeg &SegTimes, const MatNC &A0, const MatNC &H0, Coef3 &CoefPhys, Fix3 &Deriv, Real &Cost) {
    // Closed-Form Minimum-Derivative QP Solve (Paper Eq. ??)
    // TRANSPILE-FLAG F14: "Paper Eq. ??" equation number missing in the source comment; preserved verbatim. Transpiled as-is per MATLAB.
    const int D = FixLeg.size(1);
    const int r = FixLeg.size(2);
    const int K = FixLeg.size(3);
    const int M = K - 1;
    const int nc = 2 * r;
    const int nv = r * K;
    POLYQP_CAPACITY(nc <= POLYQP_MAX_NC && nv <= POLYQP_MAX_NV && M <= POLYQP_MAX_SEG && SegTimes.size() >= M);

    // Assemble global endpoint-derivative Hessian Rmat over leg knots
    // Physical action integral weight T_m^(1-2r) for squared r-th derivative over [0, T_m]
    static MatNV Rmat;
    Rmat.setZero(nv, nv);
    static Array3<Real, POLYQP_MAX_NC, POLYQP_MAX_NC, POLYQP_MAX_SEG> Hseg;
    Hseg.setConstant(nc, nc, M, Real(0)); // FIX-F07: per-segment Hessians for the cancellation-free reduced assembly
    for (int m = 1; m <= M; ++m) {
        const Real Tm = V1(SegTimes, m);
        VecNC s(nc);
        for (int q = 1; q <= r; ++q) {
            V1(s, q) = std::pow(Tm, Real(q - 1));
            V1(s, r + q) = std::pow(Tm, Real(q - 1));
        }
        const MatNC s_sT = s * s.transpose();
        const MatNC Hm = ((std::pow(Tm, Real(1 - 2 * r)) * s_sT).array() * H0.array()).matrix();
        const auto ix = R1((m - 1) * r + 1, (m - 1) * r + nc);
        Rmat(ix, ix) = Rmat(ix, ix) + Hm;
        for (int j = 1; j <= nc; ++j) {
            for (int i = 1; i <= nc; ++i) Hseg(i, j, m) = M1(Hm, i, j); // FIX-F07: per-segment Hessians for the cancellation-free reduced assembly
        }
    }

    // Solve for free derivatives dP per axis: dP* = -R_PP \ (R_FP' * dF)
    Deriv.setConstant(D, r, K, Real(0));
    Cost = 0;
    for (int d = 1; d <= D; ++d) {
        static VecNV fx;
        reshapeRowSlice(FixLeg, d, 1, K, fx);
        static Eigen::Matrix<bool, Eigen::Dynamic, 1, Eigen::ColMajor, POLYQP_MAX_NV, 1> fxFixed; // NO-NAN: ~isnan(fx)
        reshapeRowSlice(FixLegFixed, d, 1, K, fxFixed);
        const FindResult<POLYQP_MAX_NV> idxF = find<POLYQP_MAX_NV>(fxFixed.array());
        const FindResult<POLYQP_MAX_NV> idxP = find<POLYQP_MAX_NV>(!fxFixed.array());

        if (idxF.numel() < r) {
            POLYQP_ERROR(IllPosedQP, "Axis %d has %d fixed BCs; at least r = %d required to pin null space.", d, idxF.numel(), r);
        }

        static VecNV dv;
        dv = fx;
        BoolRowX<POLYQP_MAX_SEG> isShiftSeg = BoolRowX<POLYQP_MAX_SEG>::Constant(1, M, false); // FIX-F07: segments of this axis whose cost is invariant to a common shift of both positions
        static MatNV Rns;
        Rns.setZero(nv, nv); // FIX-F07: Hessian of the segments that are not shift-invariant
        if (!idxP.isempty()) {
            static MatNV R_PP;
            R_PP = Real(0.5) * (Rmat(idxP.zb(), idxP.zb()) + Rmat(idxP.zb(), idxP.zb()).transpose());
            static MatNV R_PF;
            R_PF = Rmat(idxP.zb(), idxF.zb());
            static VecNV dF;
            dF = fx(idxF.zb());

            // FIX-F07: a segment with derivatives 1..r-1 fixed at both knots and both positions free has a cost invariant to a common shift; reparametrise p_{m+1} = p_m + Delta_m (d_P = Tz * z)
            static MatNV Tz;
            Tz.setIdentity(idxP.numel(), idxP.numel());
            for (int m = 1; m <= M; ++m) {
                const int iA = (m - 1) * r + 1;
                const int iB = m * r + 1;
                if (!V1(fxFixed, iA) && !V1(fxFixed, iB) && fxFixed(R1(iA + 1, iA + r - 1)).all() && fxFixed(R1(iB + 1, iB + r - 1)).all()) { // NO-NAN: isnan(fx(i)) is !fxFixed(i)
                    const int jA = findFirst(idxP, iA);
                    const int jB = findFirst(idxP, iB);
                    Tz.row(jB - 1) = Tz.row(jA - 1);
                    M1(Tz, jB, jB) = Real(1);
                    V1(isShiftSeg, m) = true;
                }
            }
            if (isShiftSeg.any()) { // FIX-F07: reduced system Tz' * R_PP * Tz, assembled so that the common-shift direction cancels exactly
                for (int m = 1; m <= M; ++m) {
                    if (!V1(isShiftSeg, m)) {
                        for (int j = 1; j <= nc; ++j) {
                            for (int i = 1; i <= nc; ++i) M1(Rns, (m - 1) * r + i, (m - 1) * r + j) = M1(Rns, (m - 1) * r + i, (m - 1) * r + j) + Hseg(i, j, m);
                        }
                    }
                }
                static MatNV Rns_PP;
                Rns_PP = Real(0.5) * (Rns(idxP.zb(), idxP.zb()) + Rns(idxP.zb(), idxP.zb()).transpose());
                static MatNV TzT_Rns_PP;
                TzT_Rns_PP.noalias() = Tz.transpose() * Rns_PP;
                R_PP.noalias() = TzT_Rns_PP * Tz;
                static MatNV R_PP_T;
                R_PP_T = R_PP.transpose();
                R_PP = Real(0.5) * (R_PP + R_PP_T);
                static MatNV Rns_PF;
                Rns_PF = Rns(idxP.zb(), idxF.zb());
                R_PF.noalias() = Tz.transpose() * Rns_PF;
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
            static VecNV d_diag;
            d_diag = R_PP.diagonal();
            // NO-NAN: the MATLAB guard also rejects a non-finite diagonal; with real inputs it is finite, so only the numerical d_diag <= 0 test remains
            if ((d_diag.array() <= 0).any()) {
                POLYQP_ERROR(SingularRPP, "Axis %d R_PP has non-positive diagonal (min = %.2e).", d, static_cast<double>(minv(d_diag)));
            }
            static VecNV s_jac;
            s_jac = msqrt(d_diag).cwiseInverse();
            static MatNV s_jac_s_jacT;
            s_jac_s_jacT.noalias() = s_jac * s_jac.transpose();
            static MatNV R_PP_eq;
            R_PP_eq = (s_jac_s_jacT.array() * R_PP.array()).matrix();
            static MatNV R_PP_eq_T;
            R_PP_eq_T = R_PP_eq.transpose();
            R_PP_eq = Real(0.5) * (R_PP_eq + R_PP_eq_T);
            static VecNV R_PF_dF;
            R_PF_dF.noalias() = R_PF * dF;
            static VecNV rhs_eq;
            rhs_eq = (s_jac.array() * R_PF_dF.array()).matrix();

            static Eigen::PartialPivLU<MatNV> lu_rc(POLYQP_MAX_NV);
            lu_rc.compute(R_PP_eq);
            const Real rc = lu_rc.rcond();
            // TRANSPILE-FLAG F07 (fixed by FIX-F07): lateral R_PP was near-singular for unbalanced segment times (cond ~ 4e13 at (10, 0.6) s); SingularRPP could trip even in double.
            if (rc < Real(1e-13)) { // NO-NAN: the isnan(rc) part of the guard was removed (real inputs give a real rcond)
                POLYQP_ERROR(SingularRPP, "Axis %d free-derivative Hessian R_PP is singular (rcond = %.2e).", d, static_cast<double>(rc));
            }
            static Eigen::LLT<MatNV> llt(POLYQP_MAX_NV);
            llt.compute(R_PP_eq);
            static VecNV R_PP_eq_rhs_eq;
            if (llt.info() == Eigen::Success) {
                R_PP_eq_rhs_eq = llt.solve(rhs_eq);
            } else {
                static Eigen::PartialPivLU<MatNV> lu(POLYQP_MAX_NV);
                lu.compute(R_PP_eq);
                R_PP_eq_rhs_eq = lu.solve(rhs_eq);
            }
            dv(idxP.zb()) = ((-R_PP_eq_rhs_eq).array() * s_jac.array()).matrix();
            if (isShiftSeg.any()) { // FIX-F07: map the reduced solution back, d_P = Tz * z
                static VecNV z;
                z = dv(idxP.zb());
                static VecNV Tz_z;
                Tz_z.noalias() = Tz * z;
                dv(idxP.zb()) = Tz_z;
            }
        }

        assignRowSlice(Deriv, d, dv);
        if (isShiftSeg.any()) { // FIX-F07: cost without the cancelling common-shift terms (a shift-invariant segment is evaluated at p_m = 0, p_{m+1} = Delta_m)
            Cost = Cost + (dv.transpose() * Rns * dv).value();
            for (int m = 1; m <= M; ++m) {
                if (V1(isShiftSeg, m)) {
                    VecNC xm = dv(R1((m - 1) * r + 1, (m - 1) * r + nc));
                    V1(xm, r + 1) = V1(xm, r + 1) - V1(xm, 1);
                    V1(xm, 1) = Real(0);
                    MatNC Hm_m(nc, nc);
                    for (int j = 1; j <= nc; ++j) {
                        for (int i = 1; i <= nc; ++i) M1(Hm_m, i, j) = Hseg(i, j, m);
                    }
                    Cost = Cost + (xm.transpose() * Hm_m * xm).value();
                }
            }
        } else {
            Cost = Cost + (dv.transpose() * Rmat * dv).value();
        }
    }

    // Recover physical-time polynomial coefficients c_i,phys = c_i,tau / T_m^i
    const MatR D0_inv = A0(R1(1, r), R1(1, r)).diagonal().cwiseInverse().asDiagonal();
    const MatR A10 = A0(R1(r + 1, nc), R1(1, r));
    const MatR A11 = A0(R1(r + 1, nc), R1(r + 1, nc));

    CoefPhys.setConstant(nc, M, D, Real(0));
    for (int m = 1; m <= M; ++m) {
        const Real Tm = V1(SegTimes, m);
        VecNC s(nc);
        for (int q = 1; q <= r; ++q) {
            V1(s, q) = std::pow(Tm, Real(q - 1));
            V1(s, r + q) = std::pow(Tm, Real(q - 1));
        }
        VecNC scalePhys(nc);
        for (int q = 1; q <= nc; ++q) V1(scalePhys, q) = std::pow(Tm, Real(q - 1));
        for (int d = 1; d <= D; ++d) {
            VecNC d_seg;
            reshapeRowSlice(Deriv, d, m, m + 1, d_seg);
            d_seg = (s.array() * d_seg.array()).matrix();
            const VecR c_top = D0_inv * d_seg(R1(1, r));
            const VecR c_bot = Eigen::PartialPivLU<MatR>(A11).solve(d_seg(R1(r + 1, nc)) - A10 * c_top);
            VecNC c_tau(nc);
            c_tau << c_top, c_bot;
            for (int q = 1; q <= nc; ++q) CoefPhys(q, m, d) = V1(c_tau, q) / V1(scalePhys, q);
        }
    }
    return PolyQPStatus::Ok;
}

}
