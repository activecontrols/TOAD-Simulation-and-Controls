function [Time, States, Inputs] = abortFallback(X_est, Constants, vCapZ)
%% FIX-F08: Closed-Form Vertical Descent Fallback (cannot fail)
    % FIX-F08: hold the lateral position, descend rest-to-rest from z0 to 0 with the degree-9 (r = 5) profile, all derivatives zero at both ends
    MinSegTime = 0.2;
    DtSample = 0.1;
    c_v = 315 / 128;

    x0 = X_est(5);
    y0 = X_est(6);
    z0 = X_est(7);
    m_lox = X_est(14);
    m_ipa = X_est(15);

    % FIX-F08: non-finite estimates are replaced so the output stays finite (position 0, propellant full tank)
    if ~isfinite(x0)
        x0 = 0;
    end
    if ~isfinite(y0)
        y0 = 0;
    end
    if ~isfinite(z0)
        z0 = 0;
    end
    if ~isfinite(m_lox)
        m_lox = Constants.OxMass;
    end
    if ~isfinite(m_ipa)
        m_ipa = Constants.FuMass;
    end

    % FIX-F08: duration from the peak velocity c_v * |z0| / T of the rest-to-rest profile at the vertical cap
    T = max(MinSegTime, c_v * abs(z0) / vCapZ);

    % FIX-F08: same output grid as sampleOutputs
    Nc = min(3999, max(1, ceil(T / DtSample)));
    Nn = Nc + 1;
    t_vec = linspace(0, T, Nn);

    % FIX-F08: hover feedforward, upright attitude, zero rates, propellant held at the estimate
    m0 = Constants.m_dry + m_lox + m_ipa;
    T_hover = m0 * Constants.g;

    States = zeros(Nn, 15);
    for k = 1:Nn
        tau = t_vec(k) / T;
        s  = 126 * tau^5 - 420 * tau^6 + 540 * tau^7 - 315 * tau^8 + 70 * tau^9;
        sd = 630 * tau^4 - 2520 * tau^5 + 3780 * tau^6 - 2520 * tau^7 + 630 * tau^8;
        States(k, 1) = 1;
        States(k, 5) = x0;
        States(k, 6) = y0;
        States(k, 7) = z0 - z0 * s;
        States(k, 10) = -z0 * sd / T;
        States(k, 14) = m_lox;
        States(k, 15) = m_ipa;
    end

    Time = t_vec(:);
    Inputs = zeros(Nn, 4);
    for k = 1:Nn
        Inputs(k, 3) = T_hover;
    end
end
