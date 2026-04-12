function [VRE, T, f] = VRECalc(throttle, G_RMAX, kGrom, bGrom, Kg2)
    % Samples
    Num = 2 * 10^2; 
    f = logspace(-2, 4, Num);
    lowEnd = 40;                    % TADPOLE chug frequency, uncertain
    highEnd = 3000;                 % uncertain
    GrmsIN = G_RMAX * throttle;
    Broadband = GrmsIN^2 / (highEnd - lowEnd);
    n = 4;
    m = 8;
    
    % Driving PSD
    inPSD = Broadband .* abs(1 ./ (1 + (lowEnd ./ f) .^ (2 * n))) .* abs( 1 ./ (1 + (f / highEnd) .^ (2 * m)));
    mBoard = 0.1;
    
    %% Resonance calculations via 1DoF Transmissibility
    f_n = 1 / (2 * pi) * sqrt(kGrom / mBoard);
    r = f ./ f_n;
    T = sqrt((1 + (2 * bGrom * r).^2) ./ ((1 - r.^2).^2 + (2 * bGrom * r).^2));
    
    % Output PSD & GRMS
    outPSD = T.^2 .* inPSD;
    GrmsIMU = sqrt(trapz(f, outPSD));
    
    % VRE Induced bias
    VRE = Kg2 * GrmsIMU^2;
end