function [K_List, C_List] = ReadGains(trajectoryName)
    % Read Gains
    m_K = readcell(".\Guidance\Trajectories\Gains\GainMatrix" + trajectoryName);
    K_List = reshape(cell2mat(m_K), [], 12, 4);

    % Read Costs
    m_C = readcell(".\Guidance\Trajectories\Gains\CostMatrix" + trajectoryName);
    C_List = reshape(cell2mat(m_C), [], 12, 12);
end