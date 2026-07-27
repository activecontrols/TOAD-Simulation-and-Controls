function K_List = GainGenerator(X_res, U_res, t_list, constantsTOAD)
%GAINGENERATOR Generate Gains for a given trajectory
%   Detailed explanation goes here

% Hand tuning for Q for now
a_weights = ones(12,1);
a_weights = a_weights / norm(a_weights);
max_x = [0.4, 0.4, 0.1, 0.6, 0.6, 0.6, 3, 3, 3, 0.2, 0.2, 0.5];
Q = eye(12) .* a_weights ./ max_x.^2;
R = diag([12, 12, 1/2400^2, 5]);

% Augment Q with integral states
 Qi = diag([0.2, 0.2, 0.2]);
 Q = [Q zeros(12,3);
      zeros(3,12) Qi];
trajectory_gen.x = X_res;
trajectory_gen.u = [U_res'; [0, 0, (X_res(end,14)+X_res(end,15))*constantsTOAD.g, 0]]';
trajectory_gen.t = t_list;
K_List = RicattiRecursion(trajectory_gen, Q, R,constantsTOAD);



end