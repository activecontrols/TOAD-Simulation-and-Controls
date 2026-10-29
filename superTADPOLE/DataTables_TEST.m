%% Data Table Extraction Script for SKIPPER Testing
% Extracts .xlsx files first, then attaches matching .mat files from 2nd DAQ
% as supplemental measurements (downsampled to match .xlsx dimensions).

addpath('DATA');
fileListXlsx = dir('DATA\*.xlsx');
fileListMat  = dir('DATA\*.mat');

Data = struct();
KGM3_TO_LBIN3 = 3.61273e-5;
PSI_TO_PA = 6894.7;
Rho_MIN = FastFluidProps("Oxygen", 'From_P_T', 550, 130).rho;

% Build a lookup list of normalized .mat filenames for matching
matNamesNorm = strings(length(fileListMat), 1);
for m = 1:length(fileListMat)
    [~, mBase, ~] = fileparts(fileListMat(m).name);
    % Normalize spaces and hyphens to underscores, uppercase
    matNamesNorm(m) = upper(replace(replace(strtrim(mBase), ' ', '_'), '-', '_'));
end

for i = 1:length(fileListXlsx)
    currentFile = fullfile(fileListXlsx(i).folder, fileListXlsx(i).name);
    flightData = readtable(currentFile);
    
    % Determine Test Type
    if contains(fileListXlsx(i).name, "HOTFIRE", 'IgnoreCase', true)
        Data(i).TestType = "HF_TB";
    elseif contains(fileListXlsx(i).name, "CV", 'IgnoreCase', true) || contains(fileListXlsx(i).name, "CL", 'IgnoreCase', true)
        if contains(fileListXlsx(i).name, "LOX", 'IgnoreCase', true)
            Data(i).TestType = "LOX CF";
        else
            Data(i).TestType = "IPA CF";
        end
    else
        Data(i).TestType = "Undetermined";
    end

    Data(i).Time = flightData.time;
    Data(i).Phase = flightData.phase;
    Data(i).TestName = fileListXlsx(i).name;

    % Primary PT Data (.xlsx)
    Data(i).Sensor.PT_TCA_Chamber = flightData.chamber_pressure;
    Data(i).Sensor.PT_OX_ValveUp = flightData.lox_valve_upstream_pressure;
    Data(i).Sensor.PT_OX_ValveDown = flightData.lox_valve_downstream_pressure;
    Data(i).Sensor.PT_OX_Venturi = flightData.lox_venturi_differential_pressure;
    Data(i).Sensor.PT_OX_Mdot = flightData.lox_mdot;
    Data(i).Sensor.PT_IPA_ValveUp = flightData.ipa_valve_upstream_pressure;
    Data(i).Sensor.PT_IPA_ValveDown = flightData.ipa_valve_downstream_pressure;
    Data(i).Sensor.PT_IPA_Venturi = flightData.ipa_venturi_differential_pressure;
    Data(i).Sensor.PT_IPA_Mdot = flightData.ipa_mdot;

    % Primary TC Data (.xlsx)
    Data(i).Sensor.TC_OX_Venturi = flightData.lox_venturi_temperature;
    Data(i).Sensor.TC_OX_Valve = flightData.lox_valve_temperature;

    % Valve Command / Control Data (.xlsx)
    Data(i).Control.Cmd.Thrust = flightData.thrust_cmd;
    Data(i).Control.Cmd.Pos = [flightData.lox_pos_cmd, flightData.ipa_pos_cmd];
    Data(i).Control.PI.Chamber = [flightData.chamber_pressure_controller_p_component, flightData.chamber_pressure_controller_i_component];
    Data(i).Control.PI.LOX = [flightData.lox_angle_controller_p_component, flightData.lox_angle_controller_i_component];
    Data(i).Control.PI.IPA = [flightData.ipa_angle_controller_p_component, flightData.ipa_angle_controller_i_component];
    Data(i).Control.OL.Mdot = [flightData.ol_lox_mdot, flightData.ol_ipa_mdot];
    Data(i).Control.OL.Angle = [flightData.ol_lox_angle, flightData.ol_ipa_angle];

    % Valve Encoder Data (.xlsx)
    Data(i).Encoder.Pos = [flightData.lox_pos, flightData.ipa_pos];
    Data(i).Encoder.Vel = [flightData.lox_vel, flightData.ipa_vel];

    % Check for a matching supplemental .mat file from Zucrow DAQ
    [~, xlsxBase, ~] = fileparts(fileListXlsx(i).name);
    xlsxBaseNorm = upper(replace(replace(strtrim(xlsxBase), ' ', '_'), '-', '_'));
    matIdx = find(matNamesNorm == xlsxBaseNorm, 1);

    if ~isempty(matIdx)
        Data(i).HasMat = true;
        Data(i).MatFileName = fileListMat(matIdx).name;
        matFilePath = fullfile(fileListMat(matIdx).folder, fileListMat(matIdx).name);
        matStruct = load(matFilePath);

        LFMB = matStruct.allData.LFMB10k;
            raw_time      = double(LFMB.time.Value(:));
            raw_pc        = double(LFMB.pt_ch_01.Value(:));
            raw_pt_ox_up  = double(LFMB.pt10k_ox_02.Value(:));
            raw_pt_ox_inj = double(LFMB.pt_ox_04.Value(:));
            raw_pt_ox_ven = double(LFMB.pt10k_ox_01.Value(:));
            raw_pt_fu_up  = double(LFMB.pt10k_rp_04.Value(:));
            raw_pt_fu_inj = double(LFMB.pt_fu_04.Value(:));
            if isfield(LFMB, 'pv10k_rp_02_fb')
                raw_pv_fu_fb = double(LFMB.pv10k_rp_02_fb.Value(:));
            else
                raw_pv_fu_fb = nan(size(raw_time));
            end
            raw_tc_ox_ven = (double(LFMB.tc10k_ox_02.Value(:)) - 32) * 5/9 + 273.15;
            raw_tc_ox_inj = (double(LFMB.tc_ox_03.Value(:)) - 32) * 5/9 + 273.15;

        % Align time bases between 2nd DAQ (.mat) and Primary DAQ (.xlsx)
        active_idx_mat  = find(raw_pc > 50);
        active_idx_xlsx = find(Data(i).Sensor.PT_TCA_Chamber > 50);

        if ~isempty(active_idx_mat) && ~isempty(active_idx_xlsx)
            % Align using chamber pressure rise (> 50 psi) across both DAQs
            t_rise_mat  = raw_time(active_idx_mat(1));
            t_rise_xlsx = Data(i).Time(active_idx_xlsx(1));
            mat_time_aligned = (raw_time - t_rise_mat) + t_rise_xlsx;
        elseif ~isempty(active_idx_mat)
            t_start = raw_time(active_idx_mat(1)) - 1.0;
            mat_time_aligned = raw_time - t_start;
        else
            mat_time_aligned = raw_time - raw_time(1) + Data(i).Time(1);
        end

        % Ensure strictly unique time points for interpolation
        [mat_time_unique, unique_idx] = unique(mat_time_aligned, 'stable');

        % Pre-smooth high-frequency .mat data before downsampling to prevent aliasing
        dt_xlsx = mean(diff(Data(i).Time));
        dt_mat  = mean(diff(mat_time_unique));
        ds_factor = max(1, round(dt_xlsx / dt_mat));
        downsampleCh = @(v) interp1(mat_time_unique, movmedian(v(unique_idx), ds_factor), Data(i).Time, 'linear', 'extrap');

        % Populate Supplemental Sensors
        Data(i).Sensor.PT_TCA_Chamber_ZU = downsampleCh(raw_pc);
        Data(i).Sensor.PT_OX_ValveUp_ZU  = downsampleCh(raw_pt_ox_up);
        Data(i).Sensor.PT_OX_Inj         = downsampleCh(raw_pt_ox_inj); % pt_ox_04 (LOX Manifold/Injector P)
        Data(i).Sensor.PT_OX_Venturi_ZU  = downsampleCh(raw_pt_ox_ven);
        Data(i).Sensor.PT_IPA_ValveUp_ZU = downsampleCh(raw_pt_fu_up);
        Data(i).Sensor.PT_IPA_Inj        = downsampleCh(raw_pt_fu_inj); % pt_fu_04 (IPA Manifold/Injector P)
        Data(i).Sensor.TC_OX_Venturi_ZU  = downsampleCh(raw_tc_ox_ven);
        Data(i).Sensor.TC_OX_Inj         = downsampleCh(raw_tc_ox_inj); % tc_ox_03 (LOX Manifold/Injector T [K])
    else
        Data(i).HasMat = false;
        Data(i).MatFileName = "";
        
        % Initialize supplemental fields with NaNs of matching dimension
        nanVec = nan(size(Data(i).Time));
        Data(i).Sensor.PT_TCA_Chamber_ZU  = nanVec;
        Data(i).Sensor.PT_OX_ValveUp_ZU   = nanVec;
        Data(i).Sensor.PT_OX_Inj          = nanVec;
        Data(i).Sensor.PT_OX_Venturi_ZU   = nanVec;
        Data(i).Sensor.PT_IPA_ValveUp_ZU  = nanVec;
        Data(i).Sensor.PT_IPA_Inj         = nanVec;
        Data(i).Sensor.TC_OX_Venturi_ZU   = nanVec;
        Data(i).Sensor.TC_OX_Inj          = nanVec;
    end
end

%% Symmetric Filtering 
DataF = Data;

% Define filter window in seconds 
targetTimeWindow = 1; 

for i = 1:length(Data)
    % Dynamically calculate window size based on this test's sampling frequency
    dt = mean(diff(Data(i).Time));
    windowSize = max(1, round(targetTimeWindow / dt)); 

    DataF(i).Time = Data(i).Time;
    DataF(i).Phase = Data(i).Phase;

    % Primary PT Data
    DataF(i).Sensor.PT_TCA_Chamber = movmedian(Data(i).Sensor.PT_TCA_Chamber, windowSize);
    DataF(i).Sensor.PT_OX_ValveUp = movmedian(Data(i).Sensor.PT_OX_ValveUp, windowSize);
    DataF(i).Sensor.PT_OX_ValveDown = movmedian(Data(i).Sensor.PT_OX_ValveDown, windowSize);
    DataF(i).Sensor.PT_OX_Venturi = movmedian(Data(i).Sensor.PT_OX_Venturi, windowSize);
    DataF(i).Sensor.PT_OX_Mdot = movmedian(Data(i).Sensor.PT_OX_Mdot, windowSize);
    DataF(i).Sensor.PT_IPA_ValveUp = movmedian(Data(i).Sensor.PT_IPA_ValveUp, windowSize);
    DataF(i).Sensor.PT_IPA_ValveDown = movmedian(Data(i).Sensor.PT_IPA_ValveDown, windowSize);
    DataF(i).Sensor.PT_IPA_Venturi = movmedian(Data(i).Sensor.PT_IPA_Venturi, windowSize);
    DataF(i).Sensor.PT_IPA_Mdot = movmedian(Data(i).Sensor.PT_IPA_Mdot, windowSize);

    % Primary TC Data
    DataF(i).Sensor.TC_OX_Venturi = movmedian(Data(i).Sensor.TC_OX_Venturi, windowSize);
    DataF(i).Sensor.TC_OX_Valve = movmedian(Data(i).Sensor.TC_OX_Valve, windowSize);

    % Supplemental 2nd DAQ (.mat) Data
    if Data(i).HasMat
        DataF(i).Sensor.PT_TCA_Chamber_ZU  = movmedian(Data(i).Sensor.PT_TCA_Chamber_ZU, windowSize);
        DataF(i).Sensor.PT_OX_ValveUp_ZU   = movmedian(Data(i).Sensor.PT_OX_ValveUp_ZU, windowSize);
        DataF(i).Sensor.PT_OX_Inj          = movmedian(Data(i).Sensor.PT_OX_Inj, windowSize);
        DataF(i).Sensor.PT_OX_Venturi_ZU   = movmedian(Data(i).Sensor.PT_OX_Venturi_ZU, windowSize);
        DataF(i).Sensor.PT_IPA_ValveUp_ZU  = movmedian(Data(i).Sensor.PT_IPA_ValveUp_ZU, windowSize);
        DataF(i).Sensor.PT_IPA_Inj         = movmedian(Data(i).Sensor.PT_IPA_Inj, windowSize);
        DataF(i).Sensor.TC_OX_Venturi_ZU   = movmedian(Data(i).Sensor.TC_OX_Venturi_ZU, windowSize);
        DataF(i).Sensor.TC_OX_Inj          = movmedian(Data(i).Sensor.TC_OX_Inj, windowSize);
    end

    % % Recompute mdot from filtered data
    % CA_ox  = 0.066 / sqrt(1 - (0.066/0.127)^2);
    % CA_ipa = 0.062 / sqrt(1 - (0.062/0.127)^2);
    % gc = 32.174 * 12;
    % rho_ox = FastFluidProps("Oxygen", 'From_P_T', DataF(i).Sensor.PT_OX_ValveUp * PSI_TO_PA, ...
    %                         DataF(i).Sensor.TC_OX_Venturi).rho * KGM3_TO_LBIN3; 
    % rho_ipa = FastFluidProps("IPA", 'From_P_T', DataF(i).Sensor.PT_IPA_ValveUp * PSI_TO_PA, 293.15).rho * KGM3_TO_LBIN3; 
    % DataF(i).Sensor.PT_OX_Mdot = CA_ox .* sqrt(2 * rho_ox .* gc .* abs(DataF(i).Sensor.PT_OX_Venturi));
    % DataF(i).Sensor.PT_IPA_Mdot = CA_ipa .* sqrt(2 * rho_ipa .* gc .* abs(DataF(i).Sensor.PT_IPA_Venturi));

    % Valve Command / Control Data
    DataF(i).Control.Cmd.Thrust = movmedian(Data(i).Control.Cmd.Thrust, windowSize);
    DataF(i).Control.Cmd.Pos = movmedian(Data(i).Control.Cmd.Pos, windowSize);
    DataF(i).Control.PI.Chamber = movmedian(Data(i).Control.PI.Chamber, windowSize);
    DataF(i).Control.PI.LOX = movmedian(Data(i).Control.PI.LOX, windowSize);
    DataF(i).Control.PI.IPA = movmedian(Data(i).Control.PI.IPA, windowSize);
    DataF(i).Control.OL.Mdot = movmedian(Data(i).Control.OL.Mdot, windowSize);
    DataF(i).Control.OL.Angle = movmedian(Data(i).Control.OL.Angle, windowSize);

    % Valve Encoder Data
    DataF(i).Encoder.Pos = movmedian(Data(i).Encoder.Pos, windowSize);
    DataF(i).Encoder.Vel = movmedian(Data(i).Encoder.Vel, windowSize);
end


%% Analysis 
% Extract measured Cv, CdA, local injector mdot, C*, and Momentum Ratio values
Trigger = [4; 1]; % [Start offset; End offset]
TimeExclude = [8.5, 11.0; 
               15.0, 17.2];

Params = struct();
for i = 1:length(Data)

    % Extract OX Cv (across the LOX valve: uses valve upstream & downstream)
    if Data(i).TestType ~= "IPA CF"
        mdot = DataF(i).Sensor.PT_OX_Mdot;
        P_tank = DataF(i).Sensor.PT_OX_ValveUp;
        P_valve = DataF(i).Sensor.PT_OX_ValveDown;
        T_valve = DataF(i).Sensor.TC_OX_Valve;
        OX_Props = FastFluidProps("Oxygen", 'From_P_T', P_tank * PSI_TO_PA, T_valve); 
        rho = OX_Props.rho * KGM3_TO_LBIN3;
        rho_water = 0.0361; 
        Params(i).CvOX = abs(mdot ./ (rho .* rho_water .* (P_tank - P_valve)).^0.5 * 60 / 231);
    end

    % Extract FU Cv (across the IPA valve: uses valve upstream & downstream)
    if Data(i).TestType ~= "LOX CF"
        mdot = DataF(i).Sensor.PT_IPA_Mdot;
        P_tank = DataF(i).Sensor.PT_IPA_ValveUp;
        P_valve = DataF(i).Sensor.PT_IPA_ValveDown;
        IPA_Props = FastFluidProps("IPA", 'From_P_T', P_valve * PSI_TO_PA, 293.17); 
        rho = IPA_Props.rho * KGM3_TO_LBIN3;
        rho_water = 0.0361; 
        Params(i).CvIPA = abs(mdot ./ (rho .* rho_water .* (P_tank - P_valve)).^0.5 * 60 / 231);
    end

    % Extract CdA OX (uses supplemental injector/manifold P & T if available)
    if Data(i).TestType ~= "IPA CF"
        gc = 32.174 * 12;
        mdot = DataF(i).Sensor.PT_OX_Mdot;
        if Data(i).HasMat
            Pc = DataF(i).Sensor.PT_TCA_Chamber_ZU;
        else
            Pc = DataF(i).Sensor.PT_TCA_Chamber;
        end
        if Data(i).HasMat && ~all(isnan(DataF(i).Sensor.PT_OX_Inj))
            P_inj = DataF(i).Sensor.PT_OX_Inj;
        else
            P_inj = DataF(i).Sensor.PT_OX_ValveDown;
        end
        if Data(i).HasMat && ~all(isnan(DataF(i).Sensor.TC_OX_Inj))
            T_inj = DataF(i).Sensor.TC_OX_Inj;
        else
            T_inj = DataF(i).Sensor.TC_OX_Valve;
        end
        OX_Props = FastFluidProps("Oxygen", 'From_P_T', P_inj * PSI_TO_PA, T_inj); 
        rho = OX_Props.rho * KGM3_TO_LBIN3;
        Params(i).CdA_OX = abs(mdot ./ ((P_inj - Pc) .* 2 * gc .* rho).^0.5);
        Params(i).CdA_OX(rho < Rho_MIN * KGM3_TO_LBIN3) = NaN;
    end

    % Extract CdA FU (uses supplemental injector/manifold P if available)
    if Data(i).TestType ~= "LOX CF"
        gc = 32.174 * 12;
        mdot = DataF(i).Sensor.PT_IPA_Mdot;
        Pc = DataF(i).Sensor.PT_TCA_Chamber;
        if Data(i).HasMat && ~all(isnan(DataF(i).Sensor.PT_IPA_Inj))
            P_inj = DataF(i).Sensor.PT_IPA_Inj;
        else
            P_inj = DataF(i).Sensor.PT_IPA_ValveDown;
        end
        IPA_Props = FastFluidProps("IPA", 'From_P_T', P_inj * PSI_TO_PA, 293.17); 
        rho = IPA_Props.rho * KGM3_TO_LBIN3;
        Params(i).CdA_IPA = abs(mdot ./ ((P_inj - Pc) .* 2 * gc .* rho).^0.5);
    end

    %% Locally Compute Injector Mdot (Incompressible Flow from Average CdA)
    t = Data(i).Time;
    inExcludedWindow = any((t >= TimeExclude(:,1)') & (t <= TimeExclude(:,2)'), 2);
    idx_steady = (t > Trigger(1)) & (t < (t(end) - Trigger(2))) & ~inExcludedWindow;

    % Local OX Mdot across injector
    if Data(i).TestType ~= "IPA CF"
        gc = 32.174 * 12;
        Pc = DataF(i).Sensor.PT_TCA_Chamber;
        
        % Switch for available direct injector pressure
        if Data(i).HasMat && ~all(isnan(DataF(i).Sensor.PT_OX_Inj))
            P_inj_ox = DataF(i).Sensor.PT_OX_Inj;
        else
            P_inj_ox = DataF(i).Sensor.PT_OX_ValveDown;
        end
        
        % Switch for available direct injector temperature
        if Data(i).HasMat && ~all(isnan(DataF(i).Sensor.TC_OX_Inj))
            T_inj_ox = DataF(i).Sensor.TC_OX_Inj;
        else
            T_inj_ox = DataF(i).Sensor.TC_OX_Valve;
        end
        
        OX_Props = FastFluidProps("Oxygen", 'From_P_T', P_inj_ox * PSI_TO_PA, T_inj_ox);
        rho_ox = OX_Props.rho * KGM3_TO_LBIN3;
        rho_ox(rho_ox < Rho_MIN * KGM3_TO_LBIN3) = NaN;
        
        targetTimeWindow = 1; 
        dt = mean(diff(Data(i).Time));
        windowSize = max(1, round(targetTimeWindow / dt)); 
        Params(i).CdA_OX_smooth = mean(movmedian(Params(i).CdA_OX, windowSize), 'omitnan');
        dP_ox = max(0, P_inj_ox - Pc);
        Params(i).OX_Mdot_LOC = Params(i).CdA_OX_smooth .* (2 * gc .* rho_ox .* dP_ox).^0.5;

        meas_ox = DataF(i).Sensor.PT_OX_Mdot(idx_steady);
        loc_ox  = Params(i).OX_Mdot_LOC(idx_steady);
        err_ox_pct = (loc_ox - meas_ox) ./ meas_ox * 100;
        err_ox_pct(isinf(err_ox_pct)) = NaN;
        Params(i).OX_Mdot_ErrPCT = err_ox_pct;
    end

    % Local IPA Mdot across injector
    if Data(i).TestType ~= "LOX CF"
        gc = 32.174 * 12;
        Pc = DataF(i).Sensor.PT_TCA_Chamber;
        
        % Switch for available direct injector pressure
        if Data(i).HasMat && ~all(isnan(DataF(i).Sensor.PT_IPA_Inj))
            P_inj_ipa = DataF(i).Sensor.PT_IPA_Inj;
        else
            P_inj_ipa = DataF(i).Sensor.PT_IPA_ValveDown;
        end
        
        % Switch for available direct injector temperature (if field exists)
        if isfield(DataF(i).Sensor, 'TC_IPA_Inj') && ~all(isnan(DataF(i).Sensor.TC_IPA_Inj))
            T_inj_ipa = DataF(i).Sensor.TC_IPA_Inj;
        else
            T_inj_ipa = 293.17;
        end
        
        IPA_Props = FastFluidProps("IPA", 'From_P_T', P_inj_ipa * PSI_TO_PA, T_inj_ipa);
        rho_ipa = IPA_Props.rho * KGM3_TO_LBIN3;
        
        targetTimeWindow = 2; 
        dt = mean(diff(Data(i).Time));
        windowSize = max(1, round(targetTimeWindow / dt)); 
        Params(i).CdA_IPA_smooth = mean(movmedian(Params(i).CdA_IPA, windowSize), 'omitnan');
        dP_ipa = max(0, P_inj_ipa - Pc);
        Params(i).IPA_Mdot_LOC = Params(i).CdA_IPA_smooth .* (2 * gc .* rho_ipa .* dP_ipa).^0.5;

        meas_ipa = DataF(i).Sensor.PT_IPA_Mdot(idx_steady);
        loc_ipa  = Params(i).IPA_Mdot_LOC(idx_steady);
        err_ipa_pct = (loc_ipa - meas_ipa) ./ meas_ipa * 100;
        err_ipa_pct(isinf(err_ipa_pct)) = NaN;
        Params(i).IPA_Mdot_ErrPCT = (err_ipa_pct);
    end

    % Extract Cstar for hotfires (using lower-noise local injector mdot)
    if contains(DataF(i).TestType, "HF")
        g0 = 32.17;
        mdot_ox  = Params(i).OX_Mdot_LOC;
        mdot_ipa = Params(i).IPA_Mdot_LOC;
        mdot_ox = DataF(i).Sensor.PT_OX_Mdot;
        mdot_ipa = DataF(i).Sensor.PT_IPA_Mdot;
        mdot_tot = mdot_ipa + mdot_ox;
        if Data(i).HasMat
            Pc = DataF(i).Sensor.PT_TCA_Chamber_ZU;
        else
            Pc = DataF(i).Sensor.PT_TCA_Chamber;
        end
        At = 1.679;
        Params(i).Cstar = Pc .* At ./ mdot_tot .* g0;

        % Cstar Eff
        Params(i).OF    = mdot_ox ./ mdot_ipa;
        Params(i).Cstar = Pc .* At ./ mdot_tot .* g0;

        % Interpolate theoretical C* and get C* efficiency
        [Params(i).EtaCstar, Params(i).Cstar_theo] = getCstarEff(Pc, Params(i).OF, Params(i).Cstar);
    end

    % Extract Momentum Ratio for hotfires (using local injector mdot & average CdA)
    if contains(DataF(i).TestType, "HF")
        mdot_ox  = Params(i).OX_Mdot_LOC;
        mdot_ipa = Params(i).IPA_Mdot_LOC;
        mdot_ox = DataF(i).Sensor.PT_OX_Mdot;
        mdot_ipa = DataF(i).Sensor.PT_IPA_Mdot;

        % Compute densities at the injector manifold (with direct sensor switches)
        if Data(i).HasMat && ~all(isnan(DataF(i).Sensor.PT_OX_Inj))
            P_inj_ox = DataF(i).Sensor.PT_OX_Inj;
        else
            P_inj_ox = DataF(i).Sensor.PT_OX_ValveDown;
        end
        if Data(i).HasMat && ~all(isnan(DataF(i).Sensor.TC_OX_Inj))
            T_inj_ox = DataF(i).Sensor.TC_OX_Inj;
        else
            T_inj_ox = DataF(i).Sensor.TC_OX_Valve;
        end
        if Data(i).HasMat && ~all(isnan(DataF(i).Sensor.PT_IPA_Inj))
            P_inj_ipa = DataF(i).Sensor.PT_IPA_Inj;
        else
            P_inj_ipa = DataF(i).Sensor.PT_IPA_ValveDown;
        end
        if isfield(DataF(i).Sensor, 'TC_IPA_Inj') && ~all(isnan(DataF(i).Sensor.TC_IPA_Inj))
            T_inj_ipa = DataF(i).Sensor.TC_IPA_Inj;
        else
            T_inj_ipa = 293.17;
        end

        rho_ox = FastFluidProps("Oxygen", 'From_P_T', P_inj_ox * PSI_TO_PA, T_inj_ox).rho * KGM3_TO_LBIN3;
        rho_ox(rho_ox < Rho_MIN * KGM3_TO_LBIN3) = NaN;
        rho_ipa = FastFluidProps("IPA", 'From_P_T', P_inj_ipa * PSI_TO_PA, T_inj_ipa).rho * KGM3_TO_LBIN3;

        % Flow velocities (using local mdot and average injector CdA)
        Params(i).Vel_OX  = mdot_ox ./ (rho_ox .* Params(i).CdA_OX_smooth * 12);
        Params(i).Vel_IPA = mdot_ipa ./ (rho_ipa .* Params(i).CdA_IPA_smooth * 12);

        % Using effective flow areas from above
        Mom_OX  = (mdot_ox.^2) ./ (rho_ox .* Params(i).CdA_OX_smooth);
        Mom_IPA = (mdot_ipa.^2) ./ (rho_ipa .* Params(i).CdA_IPA_smooth);

        % Momentum Ratio (OX / Fuel)
        Params(i).TMR = Mom_OX ./ Mom_IPA;
        Params(i).MFR = (mdot_ox.^2 ./ rho_ox) / (mdot_ipa.^2 ./ rho_ipa);
    end

    fprintf('Test %i / %i done! (HasMat = %d)\n', i, length(Data), Data(i).HasMat);
end

%% Save Processed Workspace
save('DataTables_TEST.mat', 'Data', 'DataF', 'Params');

%% Plotting
% Generate consistent color palette and clean test labels for legends
colors = lines(length(Data));
testNames = strings(length(Data), 1);
for i = 1:length(Data)
    if Data(i).HasMat
        testNames(i) = Data(i).TestName + " (+MAT)";
    else
        testNames(i) = Data(i).TestName;
    end
end

% OX & IPA Cv vs. Valve Angle
figure('Name', 'Cv vs Angle');
for i = 1:length(Data)
    idx = (Data(i).Time > Trigger(1)) & (Data(i).Time < (Data(i).Time(end) - Trigger(2)));
    if Data(i).TestType ~= "IPA CF"
        subplot(2,1,1);
        plot(Data(i).Encoder.Pos(idx, 1) * 360, Params(i).CvOX(idx), 'x', ...
            'Color', colors(i,:), 'DisplayName', testNames(i)); 
        hold on; grid on;
    end
    if Data(i).TestType ~= "LOX CF"
        subplot(2,1,2);
        plot(Data(i).Encoder.Pos(idx, 2) * 360, Params(i).CvIPA(idx), 'x', ...
            'Color', colors(i,:), 'DisplayName', testNames(i)); 
        hold on; grid on;
    end
end
subplot(2,1,1);
xlabel('LOX Valve Angle (deg)');
ylabel('Cv_{OX}');
ylim([0, 3.2]);
xlim([20, 90]);
title('LOX Valve Cv vs. Position');
legend('show', 'Location', 'best');

subplot(2,1,2);
xlabel('IPA Valve Angle (deg)');
ylabel('Cv_{IPA}');
ylim([0, 3.2]);
xlim([20, 90]);
title('IPA Valve Cv vs. Position');
legend('show', 'Location', 'best');
sgtitle('Flow Coefficient vs. Valve Angle');

% C* vs. Chamber Pressure
figure('Name', 'Cstar vs Chamber Pressure');
for i = 1:length(Data)
    idx = (Data(i).Time > Trigger(1)) & (Data(i).Time < (Data(i).Time(end) - Trigger(2)));
    if contains(DataF(i).TestType, "HF")
        plot(DataF(i).Sensor.PT_TCA_Chamber(idx), Params(i).Cstar(idx), 'x', ...
            'Color', colors(i,:), 'DisplayName', testNames(i)); 
        hold on; grid on;
    end
end
xlabel('Chamber Pressure (psi)');
xlim([80, 300]);
ylabel('C^* (ft/s)');
title('Hotfire C^* vs. Chamber Pressure');
legend('show', 'Location', 'best');

% Injector CdA vs. Time (Solid = OX, Dashed = IPA, Color = Test)
figure('Name', 'Injector CdA vs Time');
for i = 1:length(Data)
    idx = (Data(i).Time > Trigger(1)) & (Data(i).Time < (Data(i).Time(end) - Trigger(2)));
    if Data(i).TestType ~= "IPA CF"
        plot(DataF(i).Time(idx), Params(i).CdA_OX(idx), '-', ...
            'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', testNames(i) + " (OX)"); 
        hold on; grid on;
    end
    if Data(i).TestType ~= "LOX CF"
        plot(DataF(i).Time(idx), Params(i).CdA_IPA(idx), '--', ...
            'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', testNames(i) + " (IPA)"); 
        hold on; grid on;
    end
end
ylim([0, 0.05]);
xlabel('Time (s)');
ylabel('CdA (in^2)');
title('Injector CdA vs. Time (Solid = OX, Dashed = IPA)');
legend('show', 'Location', 'best');

% Supplemental Manifold vs. Valve Downstream Comparison (For Matched Tests)
if any([Data.HasMat])
    figure('Name', 'Manifold vs Valve Downstream (Matched Tests)');
    for i = 1:length(Data)
        if Data(i).HasMat
            subplot(3,1,1);
            plot(DataF(i).Time, DataF(i).Sensor.PT_OX_ValveDown, '--', 'Color', colors(i,:), 'DisplayName', Data(i).TestName + " (Valve Down)"); hold on; grid on;
            plot(DataF(i).Time, DataF(i).Sensor.PT_OX_Inj, '-', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', Data(i).TestName + " (Manifold pt\_ox\_04)");
            
            subplot(3,1,2);
            plot(DataF(i).Time, DataF(i).Sensor.PT_IPA_ValveDown, '--', 'Color', colors(i,:), 'DisplayName', Data(i).TestName + " (Valve Down)"); hold on; grid on;
            plot(DataF(i).Time, DataF(i).Sensor.PT_IPA_Inj, '-', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', Data(i).TestName + " (Manifold pt\_fu\_04)");
            
            subplot(3,1,3);
            plot(DataF(i).Time, DataF(i).Sensor.TC_OX_Valve, '--', 'Color', colors(i,:), 'DisplayName', Data(i).TestName + " (TC Valve)"); hold on; grid on;
            plot(DataF(i).Time, DataF(i).Sensor.TC_OX_Inj, '-', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', Data(i).TestName + " (TC Manifold tc\_ox\_03)");
        end
    end
    subplot(3,1,1); ylabel('LOX Pressure (psi)'); title('LOX Valve Downstream vs. Injector Manifold Pressure'); legend('show', 'Location', 'best');
    subplot(3,1,2); ylabel('IPA Pressure (psi)'); title('IPA Valve Downstream vs. Injector Manifold Pressure'); legend('show', 'Location', 'best');
    subplot(3,1,3); xlabel('Time (s)'); ylabel('LOX Temp (K)'); title('LOX Valve vs. Injector Manifold Temperature'); legend('show', 'Location', 'best');
end

% Momentum Ratio & C* vs. Time for hotfires
figure('Name', 'Momentum Ratio vs Time');
for i = 1:length(Data)
    idx = (Data(i).Time > Trigger(1)) & (Data(i).Time < (Data(i).Time(end) - Trigger(2)));
    if contains(DataF(i).TestType, "HF")
        subplot(1,2,1);
        plot(DataF(i).Time(idx), Params(i).TMR(idx), ...
            'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', testNames(i)); 
        hold on; grid on;
        
        subplot(1,2,2);
        plot(DataF(i).Time(idx), Params(i).EtaCstar(idx), ...
            'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', testNames(i)); 
        hold on; grid on;
    end
end
subplot(1,2,1);
xlabel('Time (s)');
ylabel('Momentum Ratio (J_{OX} / J_{IPA})');
title('Hotfire Momentum Ratio vs. Time');
ylim([0, 8]);
legend('show', 'Location', 'best');

subplot(1,2,2);
xlabel('Time (s)');
ylabel('\eta C^*');
title('Hotfire \eta C^* vs. Time');
legend('show', 'Location', 'best');

% TMR vs. C* (Steady Time Gate)
figure('Name', 'TMR vs. Cstar');
for i = 1:length(Data)
    t = Data(i).Time;
    inExcludedWindow = any((t >= TimeExclude(:,1)') & (t <= TimeExclude(:,2)'), 2);
    idx = (t > Trigger(1)) & (t < (t(end) - Trigger(2))) & ~inExcludedWindow;
      
    if contains(DataF(i).TestType, "HF")
        plot(Params(i).TMR(idx), Params(i).EtaCstar(idx), 'x', ...
            'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', testNames(i)); 
        hold on; grid on;
    end
end
xlabel('Momentum Ratio (J_{OX} / J_{IPA})');
ylabel('\eta C^*');
title('Hotfire TMR vs. \eta C^* (Time-Gated)');
legend('show', 'Location', 'best');

% Combined Figure: Left = TMR vs. C*, Right = Split Velocity Scatter Colormaps
figure('Name', 'TMR and Injector Velocities vs. Cstar');

% --- LEFT COLUMN: TMR vs. C* (Spans rows 1 & 3 of a 2x2 grid) ---
subplot(2, 2, [1, 3]);
all_TMR = [];
all_Vel_IPA = [];
all_Vel_OX  = [];
all_CstarEff = [];
all_Pc = [];
all_MFR = [];

for i = 1:length(Data)
    t = Data(i).Time;
    inExcludedWindow = any((t >= TimeExclude(:,1)') & (t <= TimeExclude(:,2)'), 2);
    idx = (t > Trigger(1)) & (t < (t(end) - Trigger(2))) & ~inExcludedWindow;
      
    if contains(DataF(i).TestType, "HF")
        % Pool data for the RHS scatter plots
        all_TMR      = [all_TMR;      Params(i).TMR(idx)];
        all_MFR      = [all_MFR;      Params(i).MFR(idx)];
        all_Vel_IPA  = [all_Vel_IPA;  Params(i).Vel_IPA(idx)];
        all_Vel_OX   = [all_Vel_OX;   Params(i).Vel_OX(idx)];
        all_CstarEff = [all_CstarEff; Params(i).EtaCstar(idx)];
        all_Pc       = [all_Pc;       DataF(i).Sensor.PT_TCA_Chamber(idx)];
    end
end
valid_fit = ~isnan(all_TMR) & ~isnan(all_CstarEff);
[p, S] = polyfit(all_TMR(valid_fit), all_CstarEff(valid_fit), 1);
MomRatioSweep = linspace(1, 3, 50);
plot(MomRatioSweep, polyval(p, MomRatioSweep), 'LineWidth', 2, 'Color', 'black', 'DisplayName', 'Curve Fit'); hold on; grid on;
scatter(all_TMR, all_CstarEff, 25, all_Pc / 2.5, 'filled');
colormap(turbo);
cb1 = colorbar;
cb1.Label.String = 'Thrust Percent';
fprintf('TMR vs. CstarEff R^2: %.3f\n', S.rsquared);
xlabel('Momentum Ratio (J_{OX} / J_{IPA})');
ylabel('\eta C^*');
title('Hotfire TMR vs. C^* (Time-Gated)');
legend('show', 'Location', 'southeast');

% Clean NaNs and lock shared colorbar limits across both RHS subplots
valid = ~isnan(all_Vel_IPA) & ~isnan(all_Vel_OX) & ~isnan(all_CstarEff);
all_Vel_IPA  = all_Vel_IPA(valid);
all_Vel_OX   = all_Vel_OX(valid);
all_CstarEff = all_CstarEff(valid);

c_bounds = [min(all_CstarEff), max(all_CstarEff)];
vel_threshold = 110; % Split threshold (ft/s) between low and high regimes

% High vel
ax_high = subplot(2, 2, 2);
idx_high = all_Vel_IPA >= vel_threshold;
scatter(all_Vel_IPA(idx_high), all_Vel_OX(idx_high), 25, all_CstarEff(idx_high), 'filled');
grid on;
colormap(ax_high, turbo);
clim(ax_high, c_bounds);
cb1 = colorbar;
cb1.Label.String = '\eta C^*';
xlabel('Fuel (IPA) Injection Velocity (ft/s)');
ylabel('LOX Injection Velocity (ft/s)');
title('High Injection Velocity C^* Map');

% Low vel
ax_low = subplot(2, 2, 4);
idx_low = all_Vel_IPA < vel_threshold;
scatter(all_Vel_IPA(idx_low), all_Vel_OX(idx_low), 25, all_CstarEff(idx_low), 'filled');
grid on;
colormap(ax_low, turbo);
clim(ax_low, c_bounds);
cb2 = colorbar;
cb2.Label.String = '\eta C^*';
xlabel('Fuel (IPA) Injection Velocity (ft/s)');
ylabel('LOX Injection Velocity (ft/s)');
title('Low Injection Velocity C^* Map');

figure('Name', 'Injection Vel vs. Cstar Efficiency')
scatter(all_Vel_IPA, all_Vel_OX, 25, all_CstarEff, 'filled'); grid on;
colormap(turbo);
cb1 = colorbar;
cb1.Label.String = '\eta C*';
xlabel('Fuel (IPA) Injection Velocity (ft/s)');
ylabel('GOX Injection Velocity (ft/s)');

% Time vs. OF (Steady Time Gate)
figure('Name', 'OF Ratio Traces');
for i = 1:length(Data)
    t = Data(i).Time;
    idx = (t > Trigger(1)) & (t < (t(end) - Trigger(2)));

    if contains(DataF(i).TestType, "HF")
        plot(Data(i).Time(idx), Params(i).OF(idx), 'x', ...
            'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', testNames(i)); 
        hold on; grid on;
    end
end
xlabel('Time [s]');
ylabel('OF Ratio');
title('OF Ratio vs Time Traces');
legend('show', 'Location', 'best');

% Local Injector Mdot vs. Measured Venturi Mdot
figure('Name', 'Local vs Measured Mdot');
for i = 1:length(Data)
    inExcludedWindow = any((Data(i).Time >= TimeExclude(:,1)') & (Data(i).Time <= TimeExclude(:,2)'), 2);
    idx = (Data(i).Time > Trigger(1)) & (Data(i).Time < (Data(i).Time(end) - Trigger(2))) & ~inExcludedWindow;

    if Data(i).TestType ~= "IPA CF"
        subplot(2,1,1);
        plot(DataF(i).Time(idx), Params(i).OX_Mdot_ErrPCT, '--', ...
            'Color', colors(i,:), 'LineWidth', 1.0, 'DisplayName', testNames(i) + " (Venturi)"); 
        hold on; grid on;
    end

    if Data(i).TestType ~= "LOX CF"
        subplot(2,1,2);
        plot(DataF(i).Time(idx), Params(i).IPA_Mdot_ErrPCT, '--', ...
            'Color', colors(i,:), 'LineWidth', 1.0, 'DisplayName', testNames(i) + " (Venturi)"); 
        hold on; grid on;
    end
end
subplot(2,1,1);
xlabel('Time (s)');
ylabel('LOX Mass Flow Rate Error (%)');
title('LOX Mass Flow Rate Error Measured Venturi vs. Local Injector');
legend('show', 'Location', 'best');

subplot(2,1,2);
xlabel('Time (s)');
ylabel('IPA Mass Flow Rate Error (%)');
title('IPA Mass Flow Rate Error Measured Venturi vs. Local Injector');
legend('show', 'Location', 'best');
sgtitle('Local Injector vs. Measured Venturi Mass Flow Rates');