function Score = Acquisition(xstar, X_train, LowerBounds, UpperBounds, L_LIFE, alpha_LIFE, ThetaOpt_LIFE, L_PRESS, alpha_PRESS, ThetaOpt_PRESS, MaxDP_Scaled, BestValidScaled)
    % Reality check 
    if any(xstar < LowerBounds) || any(xstar > UpperBounds)
        Score = realmax; 
        return;
    end
    
    % Call predictors with precomputed L and alpha
    [muLIFE, stdLIFE] = PredictGP(xstar, X_train, L_LIFE, alpha_LIFE, ThetaOpt_LIFE);
    [muPRESS, stdPRESS] = PredictGP(xstar, X_train, L_PRESS, alpha_PRESS, ThetaOpt_PRESS);
    
    % Calculate expected improvement
    xi = -0.0;
    if stdLIFE > 1e-9 
        Z = (muLIFE - BestValidScaled - xi) / stdLIFE;
        EI = stdLIFE * Z * NormCDF(Z) + stdLIFE * NormPDF(Z);
    else
        EI = max(0, muLIFE - BestValidScaled - xi);
    end

    % Prob of feasability
    if stdPRESS > 1e-9
        PoF = NormCDF((MaxDP_Scaled - muPRESS) / stdPRESS);
    else
        PoF = double(muPRESS <= MaxDP_Scaled);
    end
    Score = -EI * PoF;
end