function [K_List, CostList, LA, LT] = GainGenerator(X_res, U_res, t_list, constantsTOAD)
%GAINGENERATOR Generate Gains for a given trajectory
%   Detailed explanation goes here

Q = constantsTOAD.Q_Control;
R = constantsTOAD.R_Control;

trajectory_gen.x = X_res;
trajectory_gen.u = [U_res'; [0, 0, (X_res(end,14)+X_res(end,15))*constantsTOAD.g, 0]]';
trajectory_gen.t = t_list;
[K_List, CostList, LA, LT] = RicattiRecursion(trajectory_gen, Q, R,constantsTOAD);



end