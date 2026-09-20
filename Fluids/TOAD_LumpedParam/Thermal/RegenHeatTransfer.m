function [Q_regen, DeltaP_regen, T_wall] = RegenHeatTransfer(Pc_Pa, mdot_fu, OF)
% REGENHEATTRANSFER Evaluates regenerative cooling heat transfer and channel
% pressure drop for the SKIPPER engine using precomputed 2D table lookups.
%
% Inputs:
%   Pc_Pa       - Combustion chamber pressure [Pa]
%   mdot_fu     - Fuel mass flow through coolant channels [kg/s]
%   OF          - (Optional) Current mixture ratio (default: 1.2)
%
% Outputs:
%   Q_regen      - Heat transfer rate out of chamber into fuel [W]
%   DeltaP_regen - Coolant jacket hydraulic pressure drop [Pa]
%   T_wall       - Estimated hot-gas wall temperature [K]

    persistent GI_Q GI_DP GI_Tw GI_mdot Loaded
    psi2Pa = 6894.757;
    Pa2psi = 1.0 / psi2Pa;

    if nargin < 3 || isempty(OF)
        OF = 1.2;
    end

    % Initialize Interpolants on first call
    if isempty(Loaded) || ~Loaded
        matFile = which('RegenTable_Data.mat');
        if isempty(matFile)
            matPaths = {
                fullfile(fileparts(mfilename('fullpath')), 'RegenTable_Data.mat'), ...
                fullfile(pwd, 'sandbox', 'experiments', 'Thermal', 'RegenTable_Data.mat')
            };
            for mp = 1:length(matPaths)
                if exist(matPaths{mp}, 'file')
                    matFile = matPaths{mp};
                    break;
                end
            end
        end

        if ~isempty(matFile)
            Tbl = load(matFile);
            GI_Q  = griddedInterpolant({Tbl.Pc_vec, Tbl.OF_vec}, Tbl.Q_table, 'linear', 'nearest');
            GI_DP = griddedInterpolant({Tbl.Pc_vec, Tbl.OF_vec}, Tbl.DP_table, 'linear', 'nearest');
            GI_Tw = griddedInterpolant({Tbl.Pc_vec, Tbl.OF_vec}, Tbl.Tw_table, 'linear', 'nearest');
            if isfield(Tbl, 'mdot_nom_table')
                GI_mdot = griddedInterpolant({Tbl.Pc_vec, Tbl.OF_vec}, Tbl.mdot_nom_table, 'linear', 'nearest');
            else
                GI_mdot = [];
            end
            Loaded = true;
        else
            Loaded = false;
        end
    end

    Pc_psi = max(14.7, min(350.0, Pc_Pa * Pa2psi));
    OF_clamped = max(0.8, min(2.0, OF));

    if Loaded && Pc_psi > 20.0 && mdot_fu > 0.01
        % Evaluate from pre-computed table
        Q_base  = GI_Q(Pc_psi, OF_clamped);
        DP_base = GI_DP(Pc_psi, OF_clamped);
        Tw_base = GI_Tw(Pc_psi, OF_clamped);

        % Flow rate adjustment based on physical nominal mass flow
        if ~isempty(GI_mdot)
            mdot_nom = GI_mdot(Pc_psi, OF_clamped);
        else
            mdot_nom = (Pc_psi / 250.0) * 0.5887;
        end
        flowRatio = max(0.1, min(3.0, mdot_fu / max(mdot_nom, 0.05)));

        % Heat transfer scales weakly with coolant velocity (Nu ~ Re^0.8)
        Q_regen = Q_base * (flowRatio^0.25);

        % Pressure drop scales with mdot^1.75
        DeltaP_psi = DP_base * (flowRatio^1.75);
        DeltaP_regen = DeltaP_psi * psi2Pa;

        T_wall = Tw_base;
    else
        % Engine cold or unlit
        Q_regen      = 0.0;
        DeltaP_regen = max(0.0, 1000.0 * (mdot_fu / 0.5)^1.75); % minor cold friction
        T_wall       = 293.15;
    end
end
