function Out = PolyEval(QPTraj, t, k)
%% Stateless Evaluation of k-th Derivative of Physical-Time Polynomial
    if nargin < 3
        k = 0;
    end

    t_row = t(:).';
    nPts = numel(t_row);
    M = QPTraj.Count;
    D = size(QPTraj.Coef, 3);
    nc = 2 * QPTraj.R;

    Out = zeros(D, nPts);
    if k >= nc || M == 0
        return;
    end

    % Locate active segment index for each query time (favour flight/ramp segments at knot boundaries)
    kt = QPTraj.KeyTimes(1:(M + 1));
    T_seg = QPTraj.T(1:M);
    hasCruise = isfield(QPTraj, 'IsCruise');
    idx = ones(1, nPts);
    for q = 1:nPts
        tq = t_row(q);
        m_found = 1;
        for m = 1:M
            if tq > kt(m)
                m_found = m;
            elseif tq == kt(m)
                isSpecial = QPTraj.IsHold(m) || (hasCruise && QPTraj.IsCruise(m));
                if ~isSpecial
                    m_found = m;
                end
            end
        end
        idx(q) = m_found;
    end

    % Local segment time t_loc in [0, T_m]
    t_loc = t_row - kt(idx);
    for q = 1:nPts
        if t_loc(q) < 0
            t_loc(q) = 0;
        elseif t_loc(q) > T_seg(idx(q))
            t_loc(q) = T_seg(idx(q));
        end
    end

    % Evaluate k-th derivative in physical time
    for i = k:(nc - 1)
        f = factorial(i) / factorial(i - k);
        tp = t_loc .^ (i - k);
        for d = 1:D
            c_row = reshape(QPTraj.Coef(i + 1, idx, d), 1, nPts);
            Out(d, :) = Out(d, :) + f * (c_row .* tp);
        end
    end
end
