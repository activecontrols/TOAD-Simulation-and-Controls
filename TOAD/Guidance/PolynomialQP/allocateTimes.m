function [SegTimes, Cost] = allocateTimes(FixLeg, SegTimesInit, A0, H0, MinSegTime)
%% Analytic-Gradient Time Redistribution (Paper Sec. V)
    SegTimes = SegTimesInit;
    M = numel(SegTimes);
    [~, ~, J0] = solveLeg(FixLeg, SegTimes, A0, H0);
    Cost = J0;
    if M < 2
        return;
    end

    r = size(FixLeg, 2);
    nc = 2 * r;
    Ttot = sum(SegTimes);
    a_log = log(SegTimes(:));
    MaxIter = 12;

    % Derivative order vector k_i for segment endpoint vector [0..r-1, 0..r-1]
    % Under T_m^(-1) normalized segment weighting, H_m(T_m) = T_m^(k_i + k_j - 1) * H0
    k_vec = [0:(r - 1), 0:(r - 1)].';
    K_sum = bsxfun(@plus, k_vec, k_vec.');
    ExpMat = -1 + K_sum;

    for it = 1:MaxIter
        % Solve QP at current durations to get optimal endpoint derivatives d*
        T_cur = softmaxTimes(a_log, Ttot, MinSegTime);
        [~, Deriv, J0] = solveLeg(FixLeg, T_cur, A0, H0);

        % Analytic gradient dJ*/dT_m = d*^T * (dH_m/dT_m) * d* (envelope theorem)
        grad_T = zeros(M, 1);
        for m = 1:M
            Tm = T_cur(m);
            dHm_dTm = (ExpMat .* (Tm .^ (ExpMat - 1))) .* H0;
            gm = 0;
            for d = 1:size(FixLeg, 1)
                d_seg = reshape(Deriv(d, :, m:(m + 1)), nc, 1);
                gm = gm + d_seg.' * dHm_dTm * d_seg;
            end
            grad_T(m) = gm;
        end

        % Chain rule to log-durations: dJ/da_m = T_m * (dJ/dT_m - sum_j (T_j/Ttot) * dJ/dT_j)
        mean_grad = sum((T_cur(:) / Ttot) .* grad_T);
        g = T_cur(:) .* (grad_T - mean_grad);
        g = g - mean(g);

        if max(abs(g)) < 1e-8 * max(J0, 1e-12)
            break;
        end

        gn = g / max(abs(g));
        improved = false;
        steps = [0.4, 0.2, 0.1, 0.05, 0.025];
        for s_idx = 1:numel(steps)
            an = a_log - steps(s_idx) * gn;
            Tn = softmaxTimes(an, Ttot, MinSegTime);
            [~, ~, Jn] = solveLeg(FixLeg, Tn, A0, H0);
            if Jn < J0 * (1 - 1e-6)
                a_log = an;
                J0 = Jn;
                improved = true;
                break;
            end
        end

        if ~improved
            break;
        end
    end

    SegTimes = softmaxTimes(a_log, Ttot, MinSegTime);
    Cost = J0;
end

function T = softmaxTimes(a_log, Ttot, MinSegTime)
%% Map log-durations to positive segment times summing to Ttot
    w = exp(a_log - max(a_log));
    T = Ttot * (w(:).' / sum(w));
    T = max(T, MinSegTime);
    T = Ttot * (T / sum(T));
end
