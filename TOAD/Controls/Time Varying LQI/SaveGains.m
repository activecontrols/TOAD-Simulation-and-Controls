function SaveGains(filename, constantsTOAD)
    m = readmatrix("Guidance\Trajectories\"+filename + ".csv");
    t = m(:, 1)';
    x = m(:, 2:16)';
    u = m(1:end-1, 17:20)';
    
    [K_trans_List, K_rot_List, LA, LT] = GainGenerator(x, u, t, constantsTOAD);
    
    writematrix(K_trans_List, "Guidance\Trajectories\Gains\K_trans_" + filename + ".csv");
    writematrix(K_rot_List,   "Guidance\Trajectories\Gains\K_rot_" + filename + ".csv");
    writematrix(LA,           "Guidance\Trajectories\Gains\LA_" + filename + ".csv");
    writematrix(LT,           "Guidance\Trajectories\Gains\LT_" + filename + ".csv");
end