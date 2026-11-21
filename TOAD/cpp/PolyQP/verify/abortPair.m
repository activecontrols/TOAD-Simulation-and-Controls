function [S, Sc, okS, msg] = abortPair(X_est, c, maxNodes, u_cmd, a_dist)
% Verification only (never part of the flight build).
% One abort start through MATLAB Guidance/AbortGuidance.m and polyqp_mex('AbortGuidance') (the flight node capacity).
% The caller decides when to abort, so each side plans once (Phase 6: the ABORT / tAbort / abortCalculated trigger is gone).
%   S, Sc : structs with Time, States, Inputs, abortStatus (Sc also planStatus, the C++ PolyQPStatus identifier or '')
%   okS   : structural parity (shapes and abortStatus)
% STATIC: MATLAB has no node capacity. When the C++ returns CapacityExceeded, MATLAB's plan must need more than maxNodes
% output nodes, and the expected C++ output is MATLAB's closed-form fallback on the C++ grid (F34) at abortStatus 2; S is
% replaced by that expectation.
    if nargin < 4
        u_cmd = [];
    end
    if nargin < 5
        a_dist = [];
    end
    [S.Time, S.States, S.Inputs, S.abortStatus] = AbortGuidance(X_est, c, u_cmd, a_dist);
    [Sc.Time, Sc.States, Sc.Inputs, ~, Sc.abortStatus, Sc.planStatus] = polyqp_mex('AbortGuidance', X_est, c, u_cmd, a_dist);
    msg = '';
    if strcmp(Sc.planStatus, 'PolyQP:CapacityExceeded') && S.abortStatus ~= 2
        if numel(S.Time) <= maxNodes
            okS = false;
            msg = sprintf('C++ CapacityExceeded but the MATLAB plan has %d <= %d nodes', numel(S.Time), maxNodes);
            return;
        end
        msg = sprintf('capacity fallback (MATLAB plan needs %d > %d nodes)', numel(S.Time), maxNodes);
        [S.Time, S.States, S.Inputs] = abortFallbackOnGrid(X_est, c, 1, size(Sc.States, 1));
        S.abortStatus = 2;
    end
    okS = isequal(size(S.Time), size(Sc.Time)) && isequal(size(S.States), size(Sc.States)) && isequal(size(S.Inputs), size(Sc.Inputs)) ...
        && S.abortStatus == Sc.abortStatus;
    if ~okS
        msg = sprintf('shape/abortStatus mismatch: MATLAB %d nodes, status %d; C++ %d nodes, status %d (%s)', numel(S.Time), S.abortStatus, ...
            numel(Sc.Time), Sc.abortStatus, Sc.planStatus);
    end
end
