function Score = Acquisition(xstar, X_train, y_train, LowerBounds, UpperBounds, ThetaOpt_LIFE, ThetaOpt_PRESS, MaxDP_Scaled, BestValidScaled)
    % Reality check lol
    if any(xstar < LowerBounds) || any(xstar > UpperBounds)
        Score = realmax; 
        return;
    end
    % Call predictors
    [muLIFE, stdLIFE] = PredictGP(xstar, X_train, y_train(:, 1), ThetaOpt_LIFE);
    [muPRESS, stdPRESS] = PredictGP(xstar, X_train, y_train(:, 2), ThetaOpt_PRESS);
    
    % Calculate expected improvement
    if stdLIFE > 1e-9 
        Z = (muLIFE - BestValidScaled) / stdLIFE;
        EI = stdLIFE * Z * NormCDF(Z) + stdLIFE * NormPDF(Z);
    else
        EI = max(0, muLIFE - BestValidScaled);
    end

    % Prob of feasability
    if stdPRESS > 1e-9
        PoF = NormCDF((MaxDP_Scaled - muPRESS) / stdPRESS);
    else
        % Binary check if no uncertainty
        PoF = double(muPRESS <= MaxDP_Scaled);
    end
    Score = -EI * PoF;
end