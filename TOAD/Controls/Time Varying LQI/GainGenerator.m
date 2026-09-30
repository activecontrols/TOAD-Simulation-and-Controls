function [K_trans_List, K_rot_List, LA, LT] = GainGenerator(X_res, U_res, t_list, constantsTOAD)
    % Ensure X_res is 15 x N_nodes and U_res is 4 x N_ctrl
    if size(X_res, 1) ~= 15 && size(X_res, 2) == 15
        X_res = X_res';
    end
    if size(U_res, 1) ~= 4 && size(U_res, 2) == 4
        U_res = U_res';
    end
    
    trajectory_gen.x = X_res;
    m_final = constantsTOAD.m_dry + X_res(14, end) + X_res(15, end);
    u_pad = [0; 0; m_final * constantsTOAD.g; 0];
    trajectory_gen.u = [U_res, u_pad];
    trajectory_gen.t = t_list(:)';
    [K_trans_List, K_rot_List, LA, LT] = RicattiRecursion(trajectory_gen, constantsTOAD);
end
