function GenerateRegenTable()
% GENERATEREGENTABLE Generates a pre-computed 2D lookup table for SKIPPER
% regenerative cooling heat transfer (Q_regen [W]) and coolant jacket pressure
% drop (DeltaP_regen [psi]) across the chamber operating envelope (Pc, OF).
%
% Uses the optimal SKIPPER geometry from TestCall.m:
%   Number of Channels: 80 copper channels (C101)
%   Wall Thickness:     [0.0570, 0.0510, 0.0899] in (Chamber, Throat, Exit)
%   Aspect Ratio:       [2.8886, 2.9700, 2.9950]
%   Channel Width:      [0.0350, 0.0200] in (Chamber, Throat/Exit)
%
% Integrates directly using SKRegen2_Parameterized with exact choked nozzle
% mass flows matching TOAD engine operating conditions.

    thisDir = fileparts(mfilename('fullpath'));
    skipperDir = fullfile(thisDir, '..', '..', '..', 'SKIPPER', 'SKIPPER Regen');
    if ~exist(skipperDir, 'dir')
        skipperDir = fullfile(pwd, '..', 'SKIPPER', 'SKIPPER Regen');
    end

    addpath(thisDir);
    addpath(skipperDir);
    addpath(fullfile(skipperDir, 'cea'));
    addpath(fullfile(skipperDir, 'IPA Data'));
    addpath(fullfile(skipperDir, 'Material Data'));
    addpath(fullfile(skipperDir, 'Contours'));

    fprintf('====================================================\n');
    fprintf('  Generating SKIPPER Regenerative Cooling Table     \n');
    fprintf('  Using Direct SKRegen2_Parameterized Physics Solver\n');
    fprintf('  SKIPPER Directory: %s\n', skipperDir);
    fprintf('====================================================\n');

    fprintf('Loading SKIPPER baseline engine data...\n');
    origDir = pwd;
    cd(skipperDir);
    Data = LoadData();
    cd(origDir);

    % Optimal SKIPPER Geometry (TestCall.m / REGENSizer2.m)
    NC = 80;
    WT = [0.0570, 0.0510, 0.0899];
    AR = [2.8886, 2.9700, 2.9950];
    CW = [0.0350, 0.0200];

    % Operating Grid: Pc [psi] x OF ratio
    % Encompasses throttle envelope (100 psi to 280 psi) and mixture ratio envelope (0.9 to 1.5)
    Pc_vec = [100.0, 130.0, 160.0, 190.0, 220.0, 250.0, 280.0];
    OF_vec = [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5];

    N_pc = length(Pc_vec);
    N_of = length(OF_vec);

    Q_table        = zeros(N_pc, N_of); % [W]
    DP_table       = zeros(N_pc, N_of); % [psi]
    Tw_table       = zeros(N_pc, N_of); % [K]
    mdot_nom_table = zeros(N_pc, N_of); % [kg/s]

    % Nozzle throat parameters from TOAD system (r_t = 0.800 in)
    A_throat = 0.00129717; % m^2
    Cd_throat = 0.95;
    psi2Pa = 6894.757;

    % NASA CEA c* mapping for LOX / IPA
    CEA_OF    = [0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5];
    CEA_Cstar = [1350, 1460, 1540, 1600, 1640, 1670, 1690, 1700];

    fprintf('\nEvaluating Physical Regen Performance across %dx%d Grid (%d Points):\n', ...
        N_pc, N_of, N_pc * N_of);
    fprintf('-----------------------------------------------------------------------------\n');
    fprintf('  Pc [psi] |  O/F  | Throt | mdot_fu [kg/s] | Q_tot [kW] | DeltaP [psi] | Tw_max [K]\n');
    fprintf('-----------------------------------------------------------------------------\n');

    for p = 1:N_pc
        Pc_val = Pc_vec(p);
        throt = Pc_val / 250.0;
        Pc_Pa = Pc_val * psi2Pa;

        for o = 1:N_of
            OF_val = OF_vec(o);
            cstar = interp1(CEA_OF, CEA_Cstar, OF_val, 'pchip', 'extrap');

            % Exact choked mass flow through nozzle
            mdot_tot = (Pc_Pa * A_throat * Cd_throat) / cstar;
            mdot_coolant = mdot_tot / (1.0 + OF_val);

            % Direct evaluation via parameterized tridiagonal axial finite-difference solver
            [~, PressDrop, TempArray, Qtot] = SKRegen2_Parameterized(Data, NC, WT, AR, CW, ...
                0, [], throt, OF_val, mdot_coolant);

            Tw_max = max(TempArray);

            Q_table(p, o)        = Qtot;
            DP_table(p, o)       = PressDrop;
            Tw_table(p, o)       = Tw_max;
            mdot_nom_table(p, o) = mdot_coolant;

            fprintf('   %5.1f   |  %3.1f  | %5.1f%% |     %6.4f     |   %6.2f   |     %5.2f    |   %5.1f\n', ...
                Pc_val, OF_val, throt*100, mdot_coolant, Qtot/1000, PressDrop, Tw_max);
        end
    end
    fprintf('-----------------------------------------------------------------------------\n');

    % Save pre-computed table
    outDir = fullfile(pwd, 'sandbox', 'experiments', 'Thermal');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    tableFile = fullfile(outDir, 'RegenTable_Data.mat');
    save(tableFile, 'Pc_vec', 'OF_vec', 'Q_table', 'DP_table', 'Tw_table', 'mdot_nom_table', ...
        'NC', 'WT', 'AR', 'CW');
    fprintf('\nSuccessfully saved physical Regen Table to: %s\n', tableFile);
end

