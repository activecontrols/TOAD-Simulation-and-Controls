function [eta_cstar, cstar_theo] = getCstarEff(Pc, OF, cstar_meas)

persistent F_cstar

% Load CSV and construct griddedInterpolant only on the first call
if isempty(F_cstar)
    if ~isfile('CstarSweep.csv')
        error('CstarSweep.csv not found in path. Run the CEA extraction script first.');
    end
    data = readtable('CstarSweep.csv');

    % Sort table by Pc_psia (dim 1) then OF (dim 2) to guarantee strict grid order
    data = sortrows(data, {'OF', 'Pc_psia'});

    Pc_vec = unique(data.Pc_psia); 
    OF_vec = unique(data.OF);     

    % Reshape Cstar_fps into [length(Pc_vec) x length(OF_vec)] matrix
    Cstar_grid = reshape(data.Cstar_fps, length(Pc_vec), length(OF_vec));

    % Build 2D griddedInterpolant
    F_cstar = griddedInterpolant({Pc_vec, OF_vec}, Cstar_grid, 'makima', 'nearest');
end

% Evaluate theoretical C* at each (Pc, OF) point
cstar_theo = F_cstar(Pc, OF);

% Compute C* efficiency
eta_cstar = cstar_meas ./ cstar_theo;
end