%% Rotational Linear Extended State Observer (LESO)
% Estimates body angular rates and disturbance torques, outputting an
% equivalent actuator disturbance trim U_dist = [theta, phi, roll].

function U_dist = LESO_Attitude(GND, X_est, X_trg, U_trg, L_Att, constantsTOAD, t)
    U_dist = zeros(3,1);
end