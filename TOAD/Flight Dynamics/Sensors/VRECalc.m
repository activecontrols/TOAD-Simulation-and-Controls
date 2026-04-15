function [VRE, T, f] = VRECalc(throttle, G_RMAX, kGrom, bGrom, Kg2)
    % Samples: Focus only on the relevant spectrum
    Num = 500; 
    f = logspace(log10(10), log10(2500), Num); % 10 Hz to 2500 Hz
    
    % Tuned parameters based on general launch vehicle data
    lowEnd = 50;       % Start of the plateau (typical for GEVS/SpaceX)
    highEnd = 800;     % Start of the high-frequency roll-off
    GrmsIN = G_RMAX * throttle;
    
    % Slopes: n=1 and m=1 roughly approximate a +/- 6 dB/oct power slope
    n = 1; 
    m = 2;
    
    % Driving PSD shape
    Broadband = GrmsIN^2 / (highEnd - lowEnd);
    inPSD = Broadband .* abs(1 ./ (1 + (lowEnd ./ f) .^ (2 * n))) .* abs( 1 ./ (1 + (f / highEnd) .^ (2 * m)));
    
    % Scale to exact target GRMS (since the filters alter the total area)
    current_Grms = sqrt(trapz(f, inPSD));
    scale_factor = (GrmsIN / current_Grms)^2;
    inPSD = inPSD .* scale_factor;
    
    mBoard = 0.1;
    
    %% Resonance calculations via 1DoF Transmissibility
    f_n = 1 / (2 * pi) * sqrt(kGrom / mBoard);
    r = f ./ f_n;
    T = sqrt((1 + (2 * bGrom * r).^2) ./ ((1 - r.^2).^2 + (2 * bGrom * r).^2));
    
    % Output PSD & GRMS
    outPSD = T.^2 .* inPSD;
    GrmsIMU = sqrt(trapz(log(f), outPSD .* f));
    
    % VRE Induced bias
    VRE = Kg2 * GrmsIMU^2 * pi / 180;
end