%% Digital Sequential Notch Filter implementation for ASTRAv2.
%  All IMU sensors are meant to be passed to this array always. Includes a
%  GND switch option so the filter does not run on the ground when no
%  vibrations are present. 

function OUT = DigitalNF(IN, GND, THRUST, dT)    

% Initialize past output / input memory column 
persistent X1 X2 Y1 Y2
if isempty(X1)
    X1 = ones(size(IN, 1), 2) .* IN;
    X2 = X1;
    Y1 = X1;
    Y2 = X1;
end

% Center Freq. vs Thrust Tracks, Sampling Freq, and Notch Width
TRACK = [1.3525    42.6278;
         2.6867    84.8000];
fs = 1 / dT;
width1 = 35;     % Hz

% Sequential Notch Filter
if ~GND
    %% Notch #1
        % Setup first notch at constant frequency
        f0 = TRACK(1, 1) * THRUST + TRACK(1, 2);    %Hz
        w0 = 2 * pi * f0 / fs;
        r = exp(-pi * width1 / fs);

        % Calculate Normalization Gain
        G_num = 1 - 2*r*cos(w0) + r^2;
        G_den = 2 - 2*cos(w0);
        G = G_num / G_den;

        % Filter Out
        NUM = G * (IN - 2*cos(w0)*X1(:, 1) + X1(:, 2));
        DEN = 2*r*cos(w0)*Y1(:,1) - r^2*Y1(:,2);
        OUT = NUM + DEN;
    
        % Update memory for first Notch
        X1 = [IN    X1(:, 1)];
        Y1 = [OUT   Y1(:, 1)];
        IN = OUT;
    
    %% 2nd Order Butterworth Filter
        % Set up a 2nd Order Butterworth Filter
        f0 = 20;
        C = tan(pi * f0 / fs);
        A1 = 2*(C^2 - 1) / (1 + sqrt(2)*C + C^2);
        A2 = (1 - sqrt(2)*C + C^2) / (1 + sqrt(2)*C + C^2);
        B0 = C^2 / (1 + sqrt(2)*C + C^2);
        B1 = 2 * B0;
        B2 = B0;
        OUT = B0*IN + B1*X2(:, 1) + B2*X2(:,2) - A1*Y2(:,1) - A2*Y2(:,2);

        % Update memory for second Notch
        X2 = [IN    X2(:, 1)];
        Y2 = [OUT   Y2(:, 1)];
else
    % The output is just the input (passthrough)
    OUT = IN; 
    % --- Update states sequentially for a smooth switch-on ---

    % Stage 1: History is updated with IN and OUT
    X1 = [IN    X1(:, 1)];
    Y1 = [OUT   Y1(:, 1)];

    % The input to stage 2 is the output of stage 1
    IN_2 = OUT;
    OUT_2 = IN_2; % Passthrough for stage 2

    % Stage 2: History is updated
    X2 = [IN_2  X2(:, 1)];
    Y2 = [OUT_2 Y2(:, 1)];
    OUT = OUT_2;
end