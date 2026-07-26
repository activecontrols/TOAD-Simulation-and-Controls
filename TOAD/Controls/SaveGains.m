function SaveGains(filename)
if ~exist("constantsTOAD")
    LoadTOADSim;
end


m = readmatrix("Trajectories\"+filename);

t = m(:, 1)';
x = m(:, 2:16)';
u = m(1:end-1, 17:20)';
K_List = GainGenerator(x, u, t, constantsTOAD);

K_cells = cell(size(K_List,1),1);

size(K_List,1)

 writematrix(K_List, "Trajectories\GainMatrix"+filename)

for n = 1:size(K_List,1)
    K_cells(n) = {[K_List(n, :,:)]};
    
    % writecell({[K_List(n)]}, ".\Trajectories\GainMatrix"+filename, "WriteMode", "append")
end

disp(K_cells)

% writecell(K_cells, ".\Trajectories\GainMatrix"+filename)

end