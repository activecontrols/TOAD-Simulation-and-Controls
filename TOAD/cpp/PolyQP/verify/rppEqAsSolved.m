function [Req, hasShift, ReqRaw] = rppEqAsSolved(FixLeg, SegTimes, H0)
% Verification only (never part of the flight build).
% Per-axis equilibrated free-derivative Hessian R_PP_eq exactly as solveLeg.m builds it before the rcond guard:
%   Req{d}    : the matrix solveLeg factorises (FIX-F07 reduced system when a shift-invariant segment exists)
%   ReqRaw{d} : the matrix the pre-fix solveLeg factorised (no reparametrisation)
% Empty cells for axes with no free derivative or a non-positive diagonal.
    D = size(FixLeg, 1);
    r = size(FixLeg, 2);
    K = size(FixLeg, 3);
    M = K - 1;
    nc = 2 * r;
    nv = r * K;
    Rmat = zeros(nv, nv);
    for m = 1:M
        Tm = SegTimes(m);
        s = [Tm .^ (0:(r - 1)), Tm .^ (0:(r - 1))].';
        Hm = (Tm ^ (1 - 2 * r)) * (s * s.') .* H0;
        ix = (m - 1) * r + (1:nc);
        Rmat(ix, ix) = Rmat(ix, ix) + Hm;
    end
    Req = cell(1, D);
    ReqRaw = cell(1, D);
    hasShift = false(1, D);
    for d = 1:D
        fx = reshape(FixLeg(d, :, :), nv, 1);
        idxP = find(isnan(fx));
        if isempty(idxP)
            continue;
        end
        R_PP = 0.5 * (Rmat(idxP, idxP) + Rmat(idxP, idxP).');
        ReqRaw{d} = equilibrate(R_PP);
        Tz = eye(numel(idxP));
        for m = 1:M
            iA = (m - 1) * r + 1;
            iB = m * r + 1;
            if isnan(fx(iA)) && isnan(fx(iB)) && all(~isnan(fx((iA + 1):(iA + r - 1)))) && all(~isnan(fx((iB + 1):(iB + r - 1))))
                jA = find(idxP == iA);
                jB = find(idxP == iB);
                Tz(jB, :) = Tz(jA, :);
                Tz(jB, jB) = 1;
                hasShift(d) = true;
            end
        end
        if hasShift(d)
            R_PP = Tz.' * R_PP * Tz;
            R_PP = 0.5 * (R_PP + R_PP.');
        end
        Req{d} = equilibrate(R_PP);
    end
end

function R_PP_eq = equilibrate(R_PP)
    d_diag = diag(R_PP);
    if any(d_diag <= 0) || any(~isfinite(d_diag))
        R_PP_eq = [];
        return;
    end
    s_jac = 1 ./ sqrt(d_diag);
    R_PP_eq = (s_jac * s_jac.') .* R_PP;
    R_PP_eq = 0.5 * (R_PP_eq + R_PP_eq.');
end
