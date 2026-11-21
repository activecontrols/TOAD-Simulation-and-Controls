function [Time, States, Inputs] = abortFallbackOnGrid(X_est, Constants, vCapZ, Nn)
% Verification only (never part of the flight build).
% abortFallback.m evaluated on an Nn-node grid instead of its own (the closed-form profile is exact at any node).
    T = max(0.2, 315 / 128 * abs(X_est(7)) / vCapZ);
    t_vec = linspace(0, T, Nn);
    T_hover = (Constants.m_dry + X_est(14) + X_est(15)) * Constants.g;
    States = zeros(Nn, 15);
    for k = 1:Nn
        tau = t_vec(k) / T;
        s  = 126 * tau^5 - 420 * tau^6 + 540 * tau^7 - 315 * tau^8 + 70 * tau^9;
        sd = 630 * tau^4 - 2520 * tau^5 + 3780 * tau^6 - 2520 * tau^7 + 630 * tau^8;
        States(k, [1, 5, 6, 7, 10, 14, 15]) = [1, X_est(5), X_est(6), X_est(7) - X_est(7) * s, -X_est(7) * sd / T, X_est(14), X_est(15)];
    end
    Time = t_vec(:);
    Inputs = zeros(Nn, 4);
    Inputs(:, 3) = T_hover;
end
