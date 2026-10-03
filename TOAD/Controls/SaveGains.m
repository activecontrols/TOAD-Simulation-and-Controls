function SaveGains(filename, constantsTOAD)
    m = readmatrix(fullfile("Guidance", "Trajectories", filename + ".csv"));
    t = m(:, 1)';
    x = m(:, 2:16)';
    u = m(1:end-1, 17:20)';
    
    [K_trans_List, K_rot_List, LA, LT] = GainGenerator(x, u, t, constantsTOAD);
    
    writematrix(K_trans_List, fullfile("Guidance", "Trajectories", "Gains", "K_trans_" + filename + ".csv"));
    writematrix(K_rot_List,   fullfile("Guidance", "Trajectories", "Gains", "K_rot_" + filename + ".csv"));
    writematrix(LA,           fullfile("Guidance", "Trajectories", "Gains" , "LA_" + filename + ".csv"));
    writematrix(LT,           fullfile("Guidance", "Trajectories", "Gains","LT_" + filename + ".csv"));
end