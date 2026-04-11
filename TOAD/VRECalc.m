% Samples
Num = 10^6; 
f = linspace(0.01, 10000, Num);
lowEnd = 40;        % TADPOLE chug frequency, uncertain
highEnd = 3000;     % uncertain
G_RMAX = 6;         % uncertain, anywhere from 4-20
Broadband = G_RMAX^2 / (highEnd - lowEnd);
n = 2;
m = 8;

% Driving PSD
inPSD = Broadband .* abs(1 ./ (1 + (lowEnd ./ f) .^ (2 * n))) .* abs( 1 ./ (1 + (f / highEnd) .^ (2 * m)));

% Grommet params
kGrom = 20000;      % design param
bGrom = 0.1;        % design param
mBoard = 0.1;

%% Resonance calculations via 1DoF Transmissibility
f_n = 1 / (2 * pi) * sqrt(kGrom / mBoard);
r = f ./ f_n;
T = sqrt((1 + (2 * bGrom * r).^2) ./ ((1 - r.^2).^2 + (2 * bGrom * r).^2));

% Output PSD & GRMS
outPSD = T.^2 .* inPSD;
GrmsIMU = sqrt(trapz(f, outPSD));

% VRE Induced bias
Kg2 = 0.04; %deg/s/g^2 (UNCERTAIN, 0.002-0.08)
VRE = Kg2 * GrmsIMU^2;