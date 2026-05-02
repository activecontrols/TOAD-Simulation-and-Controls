function Score = Acquisition(xstar, X_train_life, X_train_press, LowerBounds, UpperBounds, L_LIFE, alpha_LIFE, ThetaOpt_LIFE, L_PRESS, alpha_PRESS, ThetaOpt_PRESS, MaxDP_Scaled, BestValidScaled)
    % All inputs are in normalized [0,1]^D space.
    % X_train_life  -- all training points (invalid ones carry penalty)
    % X_train_press -- valid-only training points for the pressure GP (FIX 2)

    % 1. Enforce bounds
    if any(xstar < LowerBounds) || any(xstar > UpperBounds)
        Score = realmax; 
        return;
    end

    % 2. Get GP Predictions
    [muLIFE,  stdLIFE]  = PredictGP(xstar, X_train_life,  L_LIFE,  alpha_LIFE,  ThetaOpt_LIFE);
    [muPRESS, stdPRESS] = PredictGP(xstar, X_train_press, L_PRESS, alpha_PRESS, ThetaOpt_PRESS);

    % 3. Smooth uncertainties to prevent fmincon from crashing
    stdLIFE  = max(real(stdLIFE),  1e-6);
    stdPRESS = max(real(stdPRESS), 1e-6);

    % 4. Expected Improvement
    % FIX 3: xi >= 0 required. Negative xi makes EI non-zero almost
    %         everywhere, eliminating discriminative acquisition signal.
    xi = 0.01;
    Z  = (muLIFE - BestValidScaled - xi) / stdLIFE;
    EI = stdLIFE * (Z * NormCDF(Z) + NormPDF(Z));

    % 5. Probability of Feasibility
    PoF = NormCDF((MaxDP_Scaled - muPRESS) / stdPRESS);

    % 6. Final score (minimised by fmincon)
    Score = -EI * PoF;
end