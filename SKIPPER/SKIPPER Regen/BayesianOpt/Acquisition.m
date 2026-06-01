function Score = Acquisition(xstar, X_train_life, X_train_press, LowerBounds, UpperBounds, L_LIFE, alpha_LIFE, ThetaOpt_LIFE, L_PRESS, alpha_PRESS, ThetaOpt_PRESS, MaxDP_Scaled, BestValidScaled, L_TEMP, alpha_TEMP, ThetaOpt_TEMP, WeightFactor)
    % Enforce bounds
    if any(xstar < LowerBounds) || any(xstar > UpperBounds)
        Score = realmax; 
        return;
    end

    % Get GP Predictions
    [muLIFE,  stdLIFE]  = PredictGP(xstar, X_train_life,  L_LIFE,  alpha_LIFE,  ThetaOpt_LIFE);
    [muPRESS, stdPRESS] = PredictGP(xstar, X_train_press, L_PRESS, alpha_PRESS, ThetaOpt_PRESS);
    [muTEMP,  ~]        = PredictGP(xstar, X_train_press, L_TEMP,  alpha_TEMP,  ThetaOpt_TEMP);

    % Smooth uncertainties
    stdLIFE  = max(real(stdLIFE),  1e-6);
    stdPRESS = max(real(stdPRESS), 1e-6);

    % Expected Improvement
    xi = 0.005;
    Z  = (muLIFE - BestValidScaled - xi) / stdLIFE;
    EI = stdLIFE * (Z * NormCDF(Z) + NormPDF(Z));

    % Probability of Feasibility (Pressure Drop)
    PoF = NormCDF((MaxDP_Scaled - muPRESS) / stdPRESS);
    BoundedTemp = NormCDF(muTEMP);

    % Final score with Soft Penalty
    % fmincon minimizes the score. A higher predicted temperature adds to the score, penalizing the geometry.
    Score = (-EI * PoF) + (WeightFactor * BoundedTemp);
end