Num = 10^6;
f = linspace(0, 10000, Num);

Broadband = 6;
lowEnd = 50;
highEnd = 3000;
n = 2;
m = 8;

inPSD = Broadband .* abs(1 ./ (1 + (lowEnd ./ f) .^ (2 * n))) .* abs( 1 ./ (1 + (f / highEnd) .^ (2 * m)));

kGrom = 20;
bGrom = 0.05;
mBoard = 0.050;

T = sqrt((1 + (2 * bGrom * f / ((0.5/pi) * sqrt(kGrom/mBoard))).^2) / ...
    (1 - (f / ((0.5/pi) * sqrt(kGrom/mBoard))).^2).^2 + ...
    (2 * bGrom * f / ((0.5/pi) * sqrt(kGrom/mBoard))).^2);

outPSD = T.^2 .* inPSD;

GrmsIMU = sqrt(trapz(f, outPSD))

Kg2 = 0.04;
VRE = Kg2 * GrmsIMU^2