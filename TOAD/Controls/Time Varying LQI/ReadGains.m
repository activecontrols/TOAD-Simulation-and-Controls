function [K_List, C_List, LA_List, LT_List] = ReadGains(trajectoryName)
    % Read Gains
    m_K = readcell(".\Guidance\Trajectories\Gains\GainMatrix" + trajectoryName);
    K_List = reshape(cell2mat(m_K), [], 12, 4);

    % Read Costs
    m_C = readcell(".\Guidance\Trajectories\Gains\CostMatrix" + trajectoryName);
    C_List = reshape(cell2mat(m_C), [], 12, 12);

    % Read LESO Gains
    m_LA = readcell(".\Guidance\Trajectories\Gains\LA_" + trajectoryName);
    m_LT = readcell(".\Guidance\Trajectories\Gains\LT_" + trajectoryName);
    LA_List = reshape(cell2mat(m_LA), [], 9, 6);
    LT_List = reshape(cell2mat(m_LT), [], 2, 1);
end