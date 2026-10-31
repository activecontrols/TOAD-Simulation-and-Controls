function SaveGains(filename, constantsTOAD)
    filepath = fullfile('Guidance', 'Trajectories', filename);
    m = readmatrix(filepath);
    t = m(:, 1)';
    x = m(:, 2:16)';
    u = m(1:end-1, 17:20)';

    [K_trans_List, K_rot_List, LA, LT] = GainGenerator(x, u, t, constantsTOAD);

    gainsDir = fullfile('Guidance', 'Trajectories', 'Gains');
    writematrix(K_trans_List, fullfile(gainsDir, "K_trans_" + filename));
    writematrix(K_rot_List,   fullfile(gainsDir, "K_rot_" + filename));
    writematrix(LA,           fullfile(gainsDir, "LA_" + filename));
    writematrix(LT,           fullfile(gainsDir, "LT_" + filename));
end