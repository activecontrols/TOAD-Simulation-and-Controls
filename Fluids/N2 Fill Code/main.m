%% fill setup
clc;
clearvars;

% initialize tank and copv parameters
[tankV, tankP(1), tankT(1), copvV, copvP(1), copvT(1), timeStep, regSetP] = ...
    initializeModel(); 
% set the initial state for tank and copv
[tankp(1), tankm(1), tankh(1), tankU(1), copvp(1), copvm(1), copvh(1), copvU(1), copvAlumT(1), copvCarbFibT(1), copvFibGlassT(1)] = ...
    initializeState(tankV, tankP(1), tankT(1), copvV, copvP(1), copvT(1));

%% control loop
i = 1;
runLoop = 1;
fillMode(1) = 1; % 1 is on, 0 is off
ambientCool = 0; % variable to trigger ambient cooling cycle

while runLoop

    if copvAlumT(i) > 333 %if copv gets to hot, pause fill and allow cool down

        while copvT(i) > 315
            
            fillMode(i) = 0;
            mdot(i) = 0;
            regOutP(i) = regOutP(i-1);
            isChoked(i) = 0;

            % update tank and copv state variables

            [tankm(i+1), tankU(i+1), tankP(i+1), tankT(i+1), tankh(i+1), tankp(i+1), ...
                copvm(i+1), copvU(i+1), copvP(i+1), copvT(i+1), copvh(i+1), copvp(i+1)] = ...
                updateState(tankm(i), tankU(i), tankV, tankh(i), ...
                copvm(i), copvU(i), copvV, mdot(i), timeStep);

            % calculate heat flow out of the COPV
            [Q_out(i), copvAlumT(i+1), copvCarbFibT(i+1), copvFibGlassT(i+1), h_out(i+1), Ra(i+1), deltaT(i+1)] = ...
                heatFlowCOPV(fillMode(i), timeStep, copvT(i+1), copvAlumT(i), copvCarbFibT(i), copvFibGlassT(i));

            % update COPV state
            [copvU(i+1), copvP(i+1), copvT(i+1), copvh(i+1)] = updateState2(Q_out(i), copvU(i+1), copvp(i+1), copvm(i+1));

            % ambient cooldown check - if target pressure reached 
            if copvP(i+1) > regSetP

                ambientCool = 1;
                i = i + 1; % index increment - proceed to cooldown loop

            else

                i = i + 1; %increment index - proceed normally

            end      
            
        end 

    else %fill normally

        fillMode(i) = 1;

        % calculate mass flow at this time step
        [mdot(i), isChoked(i), regOutP(i), pressureEq] = modelRegulator(tankP(i), tankT(i), copvP(i), regSetP, tankP(1));

        % update tank and copv state variables

        [tankm(i+1), tankU(i+1), tankP(i+1), tankT(i+1), tankh(i+1), tankp(i+1), ...
            copvm(i+1), copvU(i+1), copvP(i+1), copvT(i+1), copvh(i+1), copvp(i+1)] = ...
            updateState(tankm(i), tankU(i), tankV, tankh(i), ...
            copvm(i), copvU(i), copvV, mdot(i), timeStep);

        % calculate heat flow out of the COPV
        [Q_out(i), copvAlumT(i+1), copvCarbFibT(i+1), copvFibGlassT(i+1), h_out(i+1), Ra(i+1), deltaT(i+1)] = ... 
             heatFlowCOPV(fillMode(i), timeStep, copvT(i+1), copvAlumT(i), copvCarbFibT(i), copvFibGlassT(i));

        % update COPV state
        [copvU(i+1), copvP(i+1), copvT(i+1), copvh(i+1)] = updateState2(Q_out(i), copvU(i+1), copvp(i+1), copvm(i+1));


        % ambient cooldown check
        if copvP(i+1) > regSetP

            ambientCool = 1;
            i = i + 1; % index increment - proceed to cooldown loop
            
        else

            i = i + 1; % increment index - proceed normally

        end

    end
    
    % ambient cooling cycle - if target pressure is reached, we cool down
    % to temperature close to ambient - pressure will drop and we can
    % continue to fill until conditions are met 
    if ambientCool == 1 % cool copv down to 10 Kelvin above ambient temperature 

        while copvT(i) > (copvT(1) + 10)

            fillMode(i) = 0;
            mdot(i) = 0;
            regOutP(i) = regOutP(i-1);
            isChoked(i) = 0;

            % update tank and copv state variables

            [tankm(i+1), tankU(i+1), tankP(i+1), tankT(i+1), tankh(i+1), tankp(i+1), ...
                copvm(i+1), copvU(i+1), copvP(i+1), copvT(i+1), copvh(i+1), copvp(i+1)] = ...
                updateState(tankm(i), tankU(i), tankV, tankh(i), ...
                copvm(i), copvU(i), copvV, mdot(i), timeStep);

            % calculate heat flow out of the COPV
            [Q_out(i), copvAlumT(i+1), copvCarbFibT(i+1), copvFibGlassT(i+1), h_out(i+1), Ra(i+1), deltaT(i+1)] = ...
                heatFlowCOPV(fillMode(i), timeStep, copvT(i+1), copvAlumT(i), copvCarbFibT(i), copvFibGlassT(i));

            % update COPV state
            [copvU(i+1), copvP(i+1), copvT(i+1), copvh(i+1)] = updateState2(Q_out(i), copvU(i+1), copvp(i+1), copvm(i+1));

            i = i + 1; % index increment

        end

        ambientCool = 0;

    end 
    
    % main loop end condition - if target pressure and ambient pressure is
    % reached or pressure is equalized
    if copvP(i-1) > (regSetP) && copvT(i-1) <= (copvT(1) + 5)

        runLoop = 0;

    end

    if pressureEq %pressure equalized

        runLoop = 0;
        disp("Pressure Equalized");

        %final cool down to ambient + 5 Kelvin due to pressure equalization
        while copvT(i) > (copvT(1) + 5) 

            fillMode(i) = 0;
            mdot(i) = 0;
            regOutP(i) = regOutP(i-1);
            isChoked(i) = 0;

            % update tank and copv state variables

            [tankm(i+1), tankU(i+1), tankP(i+1), tankT(i+1), tankh(i+1), tankp(i+1), ...
                copvm(i+1), copvU(i+1), copvP(i+1), copvT(i+1), copvh(i+1), copvp(i+1)] = ...
                updateState(tankm(i), tankU(i), tankV, tankh(i), ...
                copvm(i), copvU(i), copvV, mdot(i), timeStep);

            % calculate heat flow out of the COPV
            [Q_out(i), copvAlumT(i+1), copvCarbFibT(i+1), copvFibGlassT(i+1), h_out(i+1), Ra(i+1), deltaT(i+1)] = ... 
                heatFlowCOPV(fillMode(i), timeStep, copvT(i+1), copvAlumT(i), copvCarbFibT(i), copvFibGlassT(i));

            % update COPV state
            [copvU(i+1), copvP(i+1), copvT(i+1), copvh(i+1)] = updateState2(Q_out(i), copvU(i+1), copvp(i+1), copvm(i+1));

            i = i + 1; % index increment

        end

    end 
    
end

%% post processing

%fill time
fillTime = i * timeStep;
fprintf("Fill Time: %.2f seconds\n", fillTime);

%active vs passive filling
activeFill = (nnz(fillMode) / numel(fillMode)) * fillTime;
fprintf("Active Fill Time: %.2f seconds\n", activeFill);


%time vectors for plotting
time = (0:i-1) .* timeStep;
time1 = time;
time1(end) = [];

% tank and copv pressures
figure(1);
plot(time, tankP .* 0.000145038, time, copvP .* 0.000145038);
xlabel("Time (s)");
ylabel("Pressure (psia)");
title("Tank & COPV Pressure vs Time");
legend("Tank", "COPV");
grid on;

%mdot
figure(2);
plot(time1, mdot);
xlabel("Time (s)");
ylabel("Mass Flow Rate (kg/s)");
title("Mass Flow Rate vs Time");
grid on;

%reg pressure output
figure(3);
plot(time1, regOutP .* 0.000145038);
xlabel("Time (s)");
ylabel("Regulator Outlet Pressure (psia)");
title("Regulator Outlet Pressure vs. Time");
grid on;

%temperature
figure(4);
plot(time, tankT, time, copvT);
xlabel("Time (s)");
ylabel("Temperature (K)");
title("Tank & COPV Temperature vs Time");
legend("Tank", "COPV", "Location", "best");
grid on;

%heat loss
figure(5);
plot(time1, Q_out);
xlabel("Time (s)");
ylabel("Heat Loss (J)");
title("Heat Loss at Every Time Step");
grid on;

%wall temperatures
figure(6);
plot(time, copvAlumT, time, copvCarbFibT, time, copvFibGlassT, time, copvT);
xlabel("Time (s)");
ylabel("Temperature (K)");
title("COPV Wall Temperature vs Time");
legend("Aluminum Liner", "Carbon Fiber Wrap", "Fiberglass Coating", "Bulk Nitrogen");
grid on;

%outside convective coefficients
figure(7)
plot(time, h_out); 
xlabel("time (s)")
ylabel("h_out")
grid on;