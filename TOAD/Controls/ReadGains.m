function [K_trans_List, K_rot_List, LA_List, LT_List] = ReadGains(trajectoryName)
    gainsDir = fullfile(".", "Guidance", "Trajectories", "Gains");

    m_K_trans = readcell(fullfile(gainsDir, "K_trans_" + trajectoryName + ".csv"));
    K_trans_List = reshape(cell2mat(m_K_trans), [], 3, 6);
    m_K_rot = readcell(fullfile(gainsDir, "K_rot_" + trajectoryName + ".csv"));
    K_rot_List = reshape(cell2mat(m_K_rot), [], 3, 6);
    m_LA = readcell(fullfile(gainsDir, "LA_" + trajectoryName + ".csv"));
    m_LT = readcell(fullfile(gainsDir, "LT_" + trajectoryName + ".csv"));
    LA_List = reshape(cell2mat(m_LA), [], 6, 3);
    LT_List = reshape(cell2mat(m_LT), [], 9, 6);
end