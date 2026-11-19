function [CoefPhys, Deriv, Cost] = solveLeg(FixLeg, SegTimes, A0, H0)
%% Closed-Form Minimum-Derivative QP Solve (Paper Eq. ??)
    D = size(FixLeg, 1);
    r = size(FixLeg, 2);
    K = size(FixLeg, 3);
    M = K - 1;
    nc = 2 * r;
    nv = r * K;

    % Assemble global endpoint-derivative Hessian Rmat over leg knots
    % Physical action integral weight T_m^(1-2r) for squared r-th derivative over [0, T_m]
    Rmat = zeros(nv, nv);
    Hseg = zeros(nc, nc, M); % FIX-F07: per-segment Hessians for the cancellation-free reduced assembly
    for m = 1:M
        Tm = SegTimes(m);
        s = [Tm .^ (0:(r - 1)), Tm .^ (0:(r - 1))].';
        Hm = (Tm ^ (1 - 2 * r)) * (s * s.') .* H0;
        ix = (m - 1) * r + (1:nc);
        Rmat(ix, ix) = Rmat(ix, ix) + Hm;
        Hseg(:, :, m) = Hm; % FIX-F07: per-segment Hessians for the cancellation-free reduced assembly
    end

    % Solve for free derivatives dP per axis: dP* = -R_PP \ (R_FP' * dF)
    Deriv = zeros(D, r, K);
    Cost = 0;
    for d = 1:D
        fx = reshape(FixLeg(d, :, :), nv, 1);
        idxF = find(~isnan(fx));
        idxP = find(isnan(fx));

        if numel(idxF) < r
            error('PolyQP:IllPosedQP', ...
                'Axis %d has %d fixed BCs; at least r = %d required to pin null space.', ...
                d, numel(idxF), r);
        end

        dv = fx;
        isShiftSeg = false(1, M); % FIX-F07: segments of this axis whose cost is invariant to a common shift of both positions
        Rns = zeros(nv, nv); % FIX-F07: Hessian of the segments that are not shift-invariant
        if ~isempty(idxP)
            R_PP = 0.5 * (Rmat(idxP, idxP) + Rmat(idxP, idxP).');
            R_PF = Rmat(idxP, idxF);
            dF = fx(idxF);

            % FIX-F07: a segment with derivatives 1..r-1 fixed at both knots and both positions free has a cost invariant to a common shift; reparametrise p_{m+1} = p_m + Delta_m (d_P = Tz * z)
            Tz = eye(numel(idxP));
            for m = 1:M
                iA = (m - 1) * r + 1;
                iB = m * r + 1;
                if isnan(fx(iA)) && isnan(fx(iB)) && all(~isnan(fx((iA + 1):(iA + r - 1)))) && all(~isnan(fx((iB + 1):(iB + r - 1))))
                    jA = find(idxP == iA);
                    jB = find(idxP == iB);
                    Tz(jB, :) = Tz(jA, :);
                    Tz(jB, jB) = 1;
                    isShiftSeg(m) = true;
                end
            end
            if any(isShiftSeg) % FIX-F07: reduced system Tz' * R_PP * Tz, assembled so that the common-shift direction cancels exactly
                for m = 1:M
                    if ~isShiftSeg(m)
                        ix = (m - 1) * r + (1:nc);
                        Rns(ix, ix) = Rns(ix, ix) + Hseg(:, :, m);
                    end
                end
                R_PP = Tz.' * (0.5 * (Rns(idxP, idxP) + Rns(idxP, idxP).')) * Tz;
                R_PP = 0.5 * (R_PP + R_PP.');
                R_PF = Tz.' * Rns(idxP, idxF);
                % FIX-F07: a shift-invariant segment enters only through Delta_m (its p_{m+1} row and column, local index r + 1)
                for m = 1:M
                    if isShiftSeg(m)
                        jB = find(idxP == m * r + 1);
                        R_PP(jB, jB) = R_PP(jB, jB) + Hseg(r + 1, r + 1, m);
                        for q = 1:numel(idxF)
                            iL = idxF(q) - (m - 1) * r;
                            if iL >= 1 && iL <= nc
                                R_PF(jB, q) = R_PF(jB, q) + Hseg(r + 1, iL, m);
                            end
                        end
                    end
                end
            end

            % Jacobi diagonal equilibration: D_inv = 1 ./ sqrt(diag(R_PP))
            d_diag = diag(R_PP);
            if any(d_diag <= 0) || any(isnan(d_diag)) || any(isinf(d_diag))
                error('PolyQP:SingularRPP', ...
                    'Axis %d R_PP has non-positive diagonal (min = %.2e).', d, min(d_diag));
            end
            s_jac = 1 ./ sqrt(d_diag);
            R_PP_eq = (s_jac * s_jac.') .* R_PP;
            R_PP_eq = 0.5 * (R_PP_eq + R_PP_eq.');
            rhs_eq  = s_jac .* (R_PF * dF);

            rc = rcond(R_PP_eq);
            if isnan(rc) || rc < 1e-13
                error('PolyQP:SingularRPP', ...
                    'Axis %d free-derivative Hessian R_PP is singular (rcond = %.2e).', d, rc);
            end
            dv(idxP) = -(R_PP_eq \ rhs_eq) .* s_jac;
            if any(isShiftSeg) % FIX-F07: map the reduced solution back, d_P = Tz * z
                dv(idxP) = Tz * dv(idxP);
            end
        end

        Deriv(d, :, :) = reshape(dv, 1, r, K);
        if any(isShiftSeg) % FIX-F07: cost without the cancelling common-shift terms (a shift-invariant segment is evaluated at p_m = 0, p_{m+1} = Delta_m)
            Cost = Cost + dv.' * Rns * dv;
            for m = 1:M
                if isShiftSeg(m)
                    xm = dv((m - 1) * r + (1:nc));
                    xm(r + 1) = xm(r + 1) - xm(1);
                    xm(1) = 0;
                    Cost = Cost + xm.' * Hseg(:, :, m) * xm;
                end
            end
        else
            Cost = Cost + dv.' * Rmat * dv;
        end
    end

    % Recover physical-time polynomial coefficients c_i,phys = c_i,tau / T_m^i
    D0_inv = diag(1 ./ diag(A0(1:r, 1:r)));
    A10 = A0((r + 1):nc, 1:r);
    A11 = A0((r + 1):nc, (r + 1):nc);

    CoefPhys = zeros(nc, M, D);
    for m = 1:M
        Tm = SegTimes(m);
        s = [Tm .^ (0:(r - 1)), Tm .^ (0:(r - 1))].';
        scalePhys = (Tm .^ (0:(nc - 1))).';
        for d = 1:D
            d_seg = s .* reshape(Deriv(d, :, m:(m + 1)), nc, 1);
            c_top = D0_inv * d_seg(1:r);
            c_bot = A11 \ (d_seg((r + 1):nc) - A10 * c_top);
            c_tau = [c_top; c_bot];
            CoefPhys(:, m, d) = c_tau ./ scalePhys;
        end
    end
end
