function K_List=ReadGains(trajectoryName)

m = readcell(".\Guidance\Trajectories\"+trajectoryName);
K_List = reshape(cell2mat(m), [],15,4);
end