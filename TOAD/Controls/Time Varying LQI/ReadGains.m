function [K_trans_List, K_rot_List, LA_List, LT_List] = ReadGains(trajectoryName)

    m_K_trans = readcell(".\Guidance\Trajectories\Gains\K_trans_" + trajectoryName + ".csv");
    K_trans_List = reshape(cell2mat(m_K_trans), [], 3, 6);
    
    m_K_rot = readcell(".\Guidance\Trajectories\Gains\K_rot_" + trajectoryName + ".csv");
    K_rot_List = reshape(cell2mat(m_K_rot), [], 3, 6);
    
    m_LA = readcell(".\Guidance\Trajectories\Gains\LA_" + trajectoryName + ".csv");
    m_LT = readcell(".\Guidance\Trajectories\Gains\LT_" + trajectoryName + ".csv");
    LA_List = reshape(cell2mat(m_LA), [], 6, 3);
    LT_List = reshape(cell2mat(m_LT), [], 9, 6);
end