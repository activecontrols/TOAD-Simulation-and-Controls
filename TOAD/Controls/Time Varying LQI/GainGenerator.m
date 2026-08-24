function [K_trans_List, K_rot_List, LA, LT] = GainGenerator(X_res, U_res, t_list, constantsTOAD)
    trajectory_gen.x = X_res;
    trajectory_gen.u = [U_res'; [0, 0, (X_res(end,14)+X_res(end,15))*constantsTOAD.g, 0]]';
    trajectory_gen.t = t_list;
    [K_trans_List, K_rot_List, LA, LT] = RicattiRecursion(trajectory_gen, constantsTOAD);
end