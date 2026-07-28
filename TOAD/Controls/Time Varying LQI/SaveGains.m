function SaveGains(filename, constantsTOAD)

m = readmatrix("Guidance\Trajectories\"+filename);

t = m(:, 1)';
x = m(:, 2:16)';
u = m(1:end-1, 17:20)';
[K_List, C_List, LA, LT] = GainGenerator(x, u, t, constantsTOAD);
K_cells = cell(size(K_List,1),1);
size(K_List,1);
writematrix(K_List, "Guidance\Trajectories\Gains\GainMatrix"+filename);
writematrix(C_List, "Guidance\Trajectories\Gains\CostMatrix"+filename);
writematrix(LA, "Guidance\Trajectories\Gains\LA_"+filename);
writematrix(LT, "Guidance\Trajectories\Gains\LT_"+filename);

for n = 1:size(K_List,1)
    K_cells(n) = {[K_List(n, :,:)]};
end

end