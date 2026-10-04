function psi = PaToImperial(Pa)
    psi = Pa * 0.0001450377;
end

% rho is kg/m^3 to lb/in^3
function imprho = rhoToImperial(rho)
    imprho = rho * 13/360000;
end

% cstar is m/s to ft/s
function impcstar = cstarToImperial(cstar)
    impcstar = cstar * 3.28084;
end

% cda is m^2 to in^2
function impcda = cdaToImperial(cda)
    impcda = cda * 39.37007874^2;
end

function impcv = cvToImperial(cv)
    impcv = cv / 15850.3 * sqrt(1.450377e-4);
end
