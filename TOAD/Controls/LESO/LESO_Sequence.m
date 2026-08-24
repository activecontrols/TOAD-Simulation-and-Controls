%% Dual LESO Sequence
% Coordinates the parallel execution of the Translational and Rotational LESOs.
% Returns a consolidated disturbance vector: Dist = [a_dist; U_dist]

function Dist = LESO_Sequence(GND, X_est, X_trg, U_trg, L_Att, L_Thrust, constantsTOAD, t)
    
    % Position/Velocity LESO (Translational Disturbance Acceleration)
    a_dist = LESO_Position(GND, X_est, X_trg, U_trg, L_Thrust, constantsTOAD, t);
    
    % Attitude/Rate LESO (Rotational Actuator Disturbance Trim)
    u_dist = LESO_Attitude(GND, X_est, X_trg, U_trg, L_Att, constantsTOAD, t);
    
    % Package Consolidated Disturbance Vector (6x1)
    Dist = [a_dist; u_dist];
end