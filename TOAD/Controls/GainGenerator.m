function K_List = GainGenerator(X_res, U_res, t_list, constantsTOAD)
%GAINGENERATOR Generate Gains for a given trajectory
%   Detailed explanation goes here

% Hand tuning for Q for now
a_weights = ones(12,1);
a_weights = a_weights / norm(a_weights);
max_x = [0.23, 0.23, 0.04, 4, 4, 4, 5, 5, 5, 0.72, 0.72, 1];
Q = eye(12) .* a_weights ./ max_x.^2;
R = diag([3.5, 3.5, 1/2400^2, 5]);

% Augment Q with integral states
% Qi = diag([0.3, 0.3, 6]);
% Q = [Q zeros(6,3);
%      zeros(3,6) Qi];
trajectory_gen.x = X_res;
trajectory_gen.u = [U_res'; [0, 0, (X_res(end,14)+X_res(end,15))*constantsTOAD.g, 0]]';
trajectory_gen.t = t_list;
K_List = RicattiRecursion(trajectory_gen, Q, R,constantsTOAD);



end