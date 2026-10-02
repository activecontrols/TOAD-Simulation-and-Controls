function PSIToImperial()

end

% rho is kg/m^3 to lb/in^3
function rhoToImperial()
    
end

% cstar is m/s to ft/s
function cstarToImperial()

end

% cda is m^2 to in^2
function cdaToImperial()

end

function impcv = cvToImperial(cv)
    impcv = cv / 15850.3 * sqrt(1.450377e-4);
end
