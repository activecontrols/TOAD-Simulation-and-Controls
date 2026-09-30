classdef PolyTrajectoryQP < handle
    %% POLYTRAJECTORYQP  Stitched minimum-derivative polynomial trajectory solved as a QP.
    %
    %   Implements the closed-form constrained QP of
    %       Richter, C., Bry, A., and Roy, N. (2016), "Polynomial Trajectory Planning
    %       for Aggressive Quadrotor Flight in Dense Indoor Environments".
    %
    %   FORMULATION
    %     * K keyframes joined by M = K-1 polynomial segments of degree 2R-1
    %       (R = 4 -> septic / minimum-snap, R = 3 -> quintic / minimum-jerk).
    %     * Every segment is expressed in NORMALISED time tau = (t - t_m)/T_m in [0,1].
    %       This makes the endpoint-derivative map A0 and the cost Hessian H0 constant,
    %       independent of segment duration, and keeps the problem well conditioned.
    %     * The decision variables are the endpoint derivatives (order 0..R-1) at every
    %       keyframe. Derivatives shared by two neighbouring segments are ONE variable,
    %       so continuity up to order R-1 holds by construction (the "stitching").
    %     * Each derivative at each keyframe on each axis is either FIXED (a number) or
    %       FREE (NaN). The QP  min d'Hd  s.t.  d_F = given  has the closed form
    %               d_P* = -H_PP \ (H_PF * d_F)                    (Richter Eq. 26)
    %       Nothing else is required: no initial guess, no iterative solver.
    %
    %   INEQUALITY / CORRIDOR CONSTRAINTS
    %     Handled as in Richter Sec. VI-B: the trajectory is sampled, the point of
    %     maximum violation of each corridor is found, a keyframe is inserted there
    %     (pinned to the corridor), the segment is split, and the QP is re-solved.
    %
    %   FIX ARRAY LAYOUT
    %     Fix is D x R x K :  Fix(axis, derivativeOrder+1, keyframe)
    %     Order 0 = position, 1 = velocity, 2 = acceleration, 3 = jerk.  NaN = free.
    %
    %   Authors: PSP Active Controls (Pablo Plata, Andrew Lullo, & Antigravity)

    properties
        Fix                     % D x R x K fixed derivatives (NaN = free)
        T                       % 1 x M segment durations [s]
        Tag                     % 1 x K cell array of keyframe labels ('' = anonymous)
        MinSegTime = 0.25       % Shortest permitted segment [s]
    end

    properties (SetAccess = private)
        R = 4                   % Derivative order that is minimised (4 = snap)
        Coef                    % 2R x M x D coefficients in normalised time
        Deriv                   % D x R x K solved keyframe derivatives
        Cost = NaN              % Sum over axes of integral (d^R x / dt^R)^2 dt
    end

    methods
        %% ------------------------------------------------------------ construction
        function obj = PolyTrajectoryQP(Fix, T, varargin)
            obj.Fix = Fix;
            obj.R   = size(Fix, 2);
            K = size(Fix, 3);
            obj.T = T(:).';
            if numel(obj.T) ~= K - 1
                error('PolyTrajectoryQP:size', 'Need K-1 = %d segment times, got %d.', K - 1, numel(obj.T));
            end
            obj.Tag = repmat({''}, 1, K);
            for i = 1:2:numel(varargin)
                obj.(varargin{i}) = varargin{i + 1};
            end
            obj.solve();
        end

        %% ------------------------------------------------------------ QP solve
        function solve(obj)
            %% SOLVE  Re-solve the closed-form QP for the current Fix / T.
            [obj.Coef, obj.Deriv, obj.Cost] = PolyTrajectoryQP.solveQP(obj.Fix, obj.T);
        end

        %% ------------------------------------------------------------ evaluation
        function out = evalAt(obj, t, k)
            %% EVALAT  k-th time derivative (k = 0 pos, 1 vel, 2 acc, 3 jerk, 4 snap ...)
            %   Returns D x numel(t). Times outside [0, T_total] are clamped.
            if nargin < 3, k = 0; end
            t = t(:).';
            M = numel(obj.T); D = size(obj.Coef, 3); nc = size(obj.Coef, 1);
            cumT = [0, cumsum(obj.T)];
            idx = sum(bsxfun(@ge, t, cumT(1:M).'), 1);
            idx = min(max(idx, 1), M);
            tau = min(max((t - cumT(idx)) ./ obj.T(idx), 0), 1);
            out = zeros(D, numel(t));
            for i = k:(nc - 1)
                f  = factorial(i) / factorial(i - k);
                Ci = reshape(obj.Coef(i + 1, idx, :), numel(t), D).';
                out = out + f * Ci .* repmat(tau .^ (i - k), D, 1);
            end
            out = out ./ repmat(obj.T(idx) .^ k, D, 1);
        end

        function kt = keyTimes(obj)
            kt = [0, cumsum(obj.T)];
        end

        function s = tagTimes(obj)
            %% TAGTIMES  Struct mapping every non-empty keyframe tag to its time.
            kt = obj.keyTimes(); s = struct();
            for k = 1:numel(obj.Tag)
                if ~isempty(obj.Tag{k}), s.(obj.Tag{k}) = kt(k); end
            end
        end

        %% ------------------------------------------------------------ topology edits
        function did_insert = insertKeyframe(obj, t, pos)
            %% INSERTKEYFRAME  Split the segment containing t; pin position (NaN = free axis).
            %   Returns true if segment was successfully split, false if t violates MinSegTime bounds.
            did_insert = false;
            kt = obj.keyTimes();
            m = find(t > kt(1:end-1) & t < kt(2:end), 1);
            if isempty(m), return; end
            
            % Guard against micro-segments that corrupt Hessian conditioning
            if (t - kt(m) < obj.MinSegTime) || (kt(m+1) - t < obj.MinSegTime)
                return;
            end
            
            slab = nan(size(obj.Fix, 1), obj.R, 1);
            slab(:, 1, 1) = pos(:);
            obj.Fix = cat(3, obj.Fix(:, :, 1:m), slab, obj.Fix(:, :, m+1:end));
            obj.Tag = [obj.Tag(1:m), {''}, obj.Tag(m+1:end)];
            obj.T   = [obj.T(1:m-1), t - kt(m), kt(m+1) - t, obj.T(m+1:end)];
            did_insert = true;
        end

        function scaleTime(obj, kappa)
            %% SCALETIME  Uniformly stretch all segment durations and re-solve.
            obj.T = obj.T * kappa;
            obj.solve();
        end

        %% ------------------------------------------------------------ time allocation
        function info = optimizeTimes(obj, maxIter)
            %% OPTIMIZETIMES  Redistribute segment durations (fixed total) to minimise the QP cost.
            %   Gradient descent on log-durations with finite-difference gradients and
            %   backtracking (Richter Sec. V). Total mission time is preserved.
            if nargin < 2, maxIter = 12; end
            M = numel(obj.T); info = struct('J0', obj.Cost, 'J', obj.Cost, 'iters', 0);
            if M < 2, return; end
            Ttot = sum(obj.T);
            fcost = @(a) PolyTrajectoryQP.timeCost(obj.Fix, a, Ttot, obj.MinSegTime);
            a = log(obj.T(:));  J0 = fcost(a);  info.J0 = J0;
            for it = 1:maxIter
                g = zeros(M, 1); h = 1e-3;
                for i = 1:M
                    ap = a; ap(i) = ap(i) + h;
                    g(i) = (fcost(ap) - J0) / h;
                end
                g = g - mean(g);
                if max(abs(g)) < 1e-9 * max(J0, 1e-12), break; end
                gn = g / max(abs(g));  improved = false;
                for s = [0.4, 0.2, 0.1, 0.05, 0.025]
                    an = a - s * gn;  Jn = fcost(an);
                    if Jn < J0 * (1 - 1e-6), a = an; J0 = Jn; improved = true; break; end
                end
                info.iters = it;
                if ~improved, break; end
            end
            w = exp(a - max(a)); Tn = Ttot * w(:).' / sum(w);
            Tn = max(Tn, obj.MinSegTime);  obj.T = Ttot * Tn / sum(Tn);
            obj.solve();  info.J = obj.Cost;
        end

        %% ------------------------------------------------------------ corridors
        function info = enforceCorridors(obj, corridors, maxInsert, dtSample)
            %% ENFORCECORRIDORS  Iterative keyframe insertion at the worst violation (Richter VI-B).
            %   corridors : cell array of structs with field  fn = @(t, r, tags) -> [viol, tgt]
            %       t    : 1 x n sample times          r    : D x n sampled positions
            %       tags : struct from tagTimes()      viol : 1 x n (>0 = violated, -Inf = inactive)
            %       tgt  : D x n keyframe positions to pin if inserted (NaN = leave axis free)
            if nargin < 3 || isempty(maxInsert), maxInsert = 30; end
            if nargin < 4 || isempty(dtSample),  dtSample  = 0.05; end
            info = struct('inserted', 0, 'max_violation', 0, 'converged', true);
            for it = 1:(maxInsert + 1)
                Ttot = sum(obj.T);
                ts = linspace(0, Ttot, max(100, ceil(Ttot / dtSample) + 1));
                rs = obj.evalAt(ts, 0);
                tg = obj.tagTimes(); kt = obj.keyTimes();
                gap = min(abs(bsxfun(@minus, kt(:), ts)), [], 1);
                ins = zeros(1 + size(rs, 1), 0);  vmax = 0;
                for c = 1:numel(corridors)
                    [viol, tgt] = corridors{c}.fn(ts, rs, tg);
                    vmax = max(vmax, max(viol));
                    v2 = viol;  v2(gap < obj.MinSegTime) = -Inf;
                    [vm, i] = max(v2);
                    if vm > 0, ins = [ins, [ts(i); tgt(:, i)]]; end %#ok<AGROW>
                end
                info.max_violation = vmax;
                if isempty(ins), break; end
                if it > maxInsert, info.converged = false; break; end
                [~, ord] = sort(ins(1, :));  ins = ins(:, ord);
                ins = ins(:, [true, diff(ins(1, :)) >= obj.MinSegTime]);
                num_new = 0;
                for j = 1:size(ins, 2)
                    if obj.insertKeyframe(ins(1, j), ins(2:end, j))
                        num_new = num_new + 1;
                    end
                end
                if num_new > 0
                    obj.solve();
                    info.inserted = info.inserted + num_new;
                else
                    break;
                end
            end
        end

        %% ------------------------------------------------------------ diagnostic queries
        function d = getDuration(obj)
            d = sum(obj.T);
        end

        function len = getPathLength(obj, numSamples)
            if nargin < 2 || isempty(numSamples), numSamples = 200; end
            ts = linspace(0, sum(obj.T), numSamples);
            vs = obj.evalAt(ts, 1);
            dt = ts(2) - ts(1);
            vmag = sqrt(sum(vs.^2, 1));
            len = sum(0.5 * (vmag(1:end-1) + vmag(2:end))) * dt;
        end

        function amax = getPeakAcceleration(obj, numSamples)
            if nargin < 2 || isempty(numSamples), numSamples = 200; end
            ts = linspace(0, sum(obj.T), numSamples);
            as = obj.evalAt(ts, 2);
            amax = max(sqrt(sum(as.^2, 1)));
        end
    end

    methods (Static)
        function Fix = makeFix(D, K, R)
            %% MAKEFIX  All-free D x R x K fixed-derivative array.
            Fix = nan(D, R, K);
        end

        function [coef, deriv, cost] = solveQP(Fix, T)
            %% SOLVEQP  Closed-form constrained QP (pure function, used inside time search).
            R = size(Fix, 2); K = size(Fix, 3); D = size(Fix, 1); M = K - 1; nc = 2 * R;
            nv = R * K;
            [A0inv, H0] = PolyTrajectoryQP.basis(R);

            % Assemble global Hessian over keyframe-major derivative vector
            % d = [d(kf1, orders 0..R-1); d(kf2, ...); ...] : segment m owns block (m-1)R+(1:2R)
            H = zeros(nv);
            for m = 1:M
                s  = repmat(T(m) .^ (0:R-1), 1, 2).';
                Hm = T(m)^(1 - 2*R) * (s * s.') .* H0;
                ix = (m-1)*R + (1:nc);
                H(ix, ix) = H(ix, ix) + Hm;
            end

            deriv = zeros(D, R, K); cost = 0;
            for d = 1:D
                fx = reshape(Fix(d, :, :), nv, 1);
                F = find(~isnan(fx)); P = find(isnan(fx));
                if numel(F) < R
                    error('PolyTrajectoryQP:underdetermined', ...
                        'Axis %d has %d fixed values; at least R = %d are required.', d, numel(F), R);
                end
                dv = fx;
                if ~isempty(P)
                    HPP = H(P, P);
                    HPP = 0.5 * (HPP + HPP.');
                    rhs = H(P, F) * fx(F);
                    rc = rcond(HPP);
                    if isnan(rc) || rc < 1e-12
                        dv(P) = -pinv(HPP) * rhs;
                    else
                        dv(P) = -(HPP \ rhs);
                    end
                end
                deriv(d, :, :) = reshape(dv, 1, R, K);
                cost = cost + dv.' * H * dv;
            end

            coef = zeros(nc, M, D);
            for m = 1:M
                s = repmat(T(m) .^ (0:R-1), 1, 2).';
                for d = 1:D
                    dseg = reshape(deriv(d, :, m:m+1), nc, 1);
                    coef(:, m, d) = A0inv * (s .* dseg);
                end
            end
        end

        function J = timeCost(Fix, a, Ttot, Tmin)
            %% TIMECOST  QP cost for log-duration vector a (softmax keeps sum(T) = Ttot).
            w = exp(a - max(a)); T = Ttot * w(:).' / sum(w);
            T = max(T, Tmin);    T = Ttot * T / sum(T);
            try
                [~, ~, J] = PolyTrajectoryQP.solveQP(Fix, T);
            catch
                J = 1e12;
            end
        end

        function [A0inv, H0] = basis(R)
            %% BASIS  Constant normalised-time matrices for polynomial degree 2R-1.
            %   A0 : coefficients -> [p^(k)(0); p^(k)(1)], k = 0..R-1   (derivatives in tau)
            %   H0 : cost Hessian in the same endpoint-derivative coordinates,
            %        H0 = A0^-T Q0 A0^-1,  Q0_ij = int_0^1 p_i^(R) p_j^(R) dtau
            persistent cR cA cH
            if ~isempty(cR) && cR == R, A0inv = cA; H0 = cH; return; end
            nc = 2 * R; A0 = zeros(nc); Q0 = zeros(nc);
            for k = 0:R-1
                A0(k + 1, k + 1) = factorial(k);
                for i = k:(nc - 1)
                    A0(R + k + 1, i + 1) = factorial(i) / factorial(i - k);
                end
            end
            for i = R:(nc - 1)
                for j = R:(nc - 1)
                    Q0(i + 1, j + 1) = (factorial(i) / factorial(i - R)) * ...
                                       (factorial(j) / factorial(j - R)) / (i + j - 2*R + 1);
                end
            end
            A0inv = inv(A0);
            H0 = A0inv.' * Q0 * A0inv;  H0 = 0.5 * (H0 + H0.');
            cR = R; cA = A0inv; cH = H0;
        end
    end
end
