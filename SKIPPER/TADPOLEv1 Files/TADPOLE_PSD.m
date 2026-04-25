% Hotfire Detection, Pressure PSD, and Estimated Vibration PSD Script
clear; clc; close all;

% 1. Set figures to dock by default
set(0, 'DefaultFigureWindowStyle', 'docked');

% --- Configuration Parameters ---
dataFolder = 'TADPOLEv1 Files';          % Set to 'Raw Data Logs' if script is outside the folder
filePattern = fullfile(dataFolder, 'TEST_*_data.mat');
testFiles = dir(filePattern);

if isempty(testFiles)
    error('No test data files found. Please check your directory path.');
end

% Detection & Clipping Thresholds
pc_min_operating = 100;   % [psi] Strict threshold: Must hit this to be a hotfire
envelope_thresh  = 50;    % [psi] Boundary threshold: Defines the start/end to survive throttles
pad_time         = 8;     % [sec] Time buffer for the VISUAL plot only

% --- Vibration Estimation Parameters ---
% Calculate the pressure-to-force coupling factor (K_pf) based on hotfire data:
thrust_cal_1 = 590 / 280; % [lbf/psi]
thrust_cal_2 = 236 / 105; % [lbf/psi]
K_pf = mean([thrust_cal_1, thrust_cal_2]); % ~2.18 lbf/psi

% Engine Mass (Update this to the actual wet mass of the suspended engine test article)
engine_mass_lbm = 12;

% GRMS Integration Bounds
grms_f_min = 10; % [Hz] Lower bound to exclude bulk fluid/throttling noise

% Data storage for the final aggregate plot
psd_data = struct('name', {}, 'f', {}, 'pxx', {}, 'pxxc', {}, 'pxx_vib', {}, 'grms', {});

for i = 1:length(testFiles)
    fileName = testFiles(i).name;
    filePath = fullfile(testFiles(i).folder, fileName);
    
    % Load Data safely
    fprintf('Evaluating %s... ', fileName);
    try
        data = load(filePath);
    catch
        fprintf('Failed to load. Skipping.\n');
        continue;
    end
    
    if ~isfield(data, 'LFMB10k')
         fprintf('Invalid data structure. Skipping.\n');
         continue;
    end
    
    t = data.LFMB10k.time.Value;
    PC = data.LFMB10k.pt_ch_01.Value;
    
    % Hotfire Detection: Did it ever reach operating pressure?
    if max(PC) < pc_min_operating
        fprintf('Cold flow / Abort. Skipping.\n');
        continue;
    end
    
    Fs = 1 / mean(diff(t));
    
    % 2. Improved Logic: Envelope Detection
    envelope_indices = find(PC > envelope_thresh);
    if isempty(envelope_indices)
        fprintf('No sustained envelope. Skipping.\n');
        continue; 
    end
    
    best_start = envelope_indices(1);
    best_end = envelope_indices(end);
    
    % Extract strict steady-state region for PSD
    trim_idx = round(0.25 * Fs);
    idx_psd_start = min(best_start + trim_idx, best_end);
    idx_psd_end   = max(best_end - trim_idx, best_start);
    
    PC_burn = PC(idx_psd_start:idx_psd_end);
    PC_zero_mean = PC_burn - mean(PC_burn);
    
    % 3. Calculate Pressure PSD with Confidence Bounds (95%)
    window = round(Fs * 0.25);
    noverlap = round(window * 0.5);
    [pxx, f, pxxc] = pwelch(PC_zero_mean, window, noverlap, [], Fs, 'ConfidenceLevel', 0.95);
    
    % --- 4. CONVERT TO VIBRATION PSD ---
    % Step A: Convert Pressure PSD [psi^2/Hz] to Force PSD [lbf^2/Hz]
    pxx_force = pxx .* (K_pf^2); 
    
    % Step B: Define Structural Transfer Function (FRF) magnitude squared [g^2/lbf^2]
    % Since resonances (>2kHz) are out of band, we assume a pure rigid body response (a = F/m)
    rigid_body_accel = 1 / engine_mass_lbm; % [g/lbf]
    FRF_mag_squared = ones(size(f)) * (rigid_body_accel^2);
    
    % Step C: Calculate Estimated Vibration PSD [g^2/Hz]
    pxx_vib = pxx_force .* FRF_mag_squared;
    
    % Step D: Calculate GRMS (filtering out low-frequency noise)
    idx_grms = find(f >= grms_f_min); % Find indices where f is >= 10 Hz
    grms_val = sqrt(trapz(f(idx_grms), pxx_vib(idx_grms)));
    fprintf('Hotfire Detected! GRMS (>%dHz) = %.3f g\n', grms_f_min, grms_val);
    
    % Store for final plot
    psd_data(end+1).name = fileName;
    psd_data(end).f = f;
    psd_data(end).pxx = pxx;
    psd_data(end).pxxc = pxxc;
    psd_data(end).pxx_vib = pxx_vib;
    psd_data(end).grms = grms_val;
    
    % 5. Visual Plotting Extraction
    pad_idx = round(pad_time * Fs);
    idx_plot_start = max(1, envelope_indices(1) - pad_idx);
    idx_plot_end = min(length(t), envelope_indices(end) + pad_idx);
    
    t_clip  = t(idx_plot_start:idx_plot_end);
    PC_clip = PC(idx_plot_start:idx_plot_end);
    
    % 6. Individual Figure Generation
    figure('Name', sprintf('Analysis: %s', fileName));
    
    % Subplot A: Clipped Autosequence
    subplot(1, 3, 1);
    plot(t_clip, PC_clip, 'r', 'LineWidth', 1.2);
    grid on; hold on;
    title(sprintf('%s - Autosequence', strrep(fileName, '_', '\_')));
    xlabel('Time [sec]'); ylabel('Pressure [psi]');
    xline(t(idx_psd_start), '--k'); xline(t(idx_psd_end), '--k');
    xlim([min(t_clip) max(t_clip)]);
    
    % Subplot B: Pressure PSD (dB/Hz)
    subplot(1, 3, 2);
    f_fill = [f; flipud(f)];
    p_fill = 10*log10([pxxc(:,1); flipud(pxxc(:,2))]);
    fill(f_fill, p_fill, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
    plot(f, 10*log10(pxx), 'b', 'LineWidth', 1.2); 
    grid on;
    title('Chamber Pressure PSD');
    xlabel('Frequency [Hz]'); ylabel('Power/Freq [dB/Hz relative to psi^2/Hz]');
    xlim([0, Fs/2]);
    
    % Subplot C: Estimated Vibration PSD (g^2/Hz)
    subplot(1, 3, 3);
    semilogy(f, pxx_vib, 'k', 'LineWidth', 1.2);
    grid on;
    title(sprintf('Estimated Vibration PSD (G_{RMS} = %.3f)', grms_val));
    xlabel('Frequency [Hz]'); ylabel('Vibration [g^2/Hz]');
    xlim([0, Fs/2]);
end

fprintf('\nBatch processing complete. Generating overlay plots...\n');

% 7. Final Aggregate Overlay Figures
if ~isempty(psd_data)
    % Figure A: Original Pressure PSD Overlay
    figure('Name', 'Aggregate Pressure PSD Overlay');
    hold on; grid on;
    color_map = lines(length(psd_data));
    plot_handles = gobjects(length(psd_data), 1);
    legend_names_psd = cell(length(psd_data), 1);
    
    for k = 1:length(psd_data)
        f_k = psd_data(k).f;
        pxx_k = 10*log10(psd_data(k).pxx);
        pxxc_k = 10*log10(psd_data(k).pxxc);
        
        f_fill = [f_k; flipud(f_k)];
        p_fill = [pxxc_k(:,1); flipud(pxxc_k(:,2))];
        fill(f_fill, p_fill, color_map(k,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
        plot_handles(k) = plot(f_k, pxx_k, 'Color', color_map(k,:), 'LineWidth', 1.5);
        legend_names_psd{k} = strrep(psd_data(k).name, '_', '\_');
    end
    title('Aggregate Chamber Pressure PSD Overlay');
    xlabel('Frequency [Hz]'); ylabel('Power/Freq [dB/Hz relative to psi^2/Hz]');
    legend(plot_handles, legend_names_psd, 'Location', 'best');
    xlim([0, psd_data(1).f(end)]);
    hold off;
    
    % Figure B: Estimated Vibration PSD Overlay with GRMS
    figure('Name', 'Aggregate Vibration PSD Overlay');
    hold on; grid on; set(gca, 'YScale', 'log');
    plot_handles_vib = gobjects(length(psd_data), 1);
    legend_names_vib = cell(length(psd_data), 1);
    
    for k = 1:length(psd_data)
        plot_handles_vib(k) = semilogy(psd_data(k).f, psd_data(k).pxx_vib, 'Color', color_map(k,:), 'LineWidth', 1.5);
        % Format legend to include the GRMS calculation
        legend_names_vib{k} = sprintf('%s (G_{RMS} = %.3f)', strrep(psd_data(k).name, '_', '\_'), psd_data(k).grms);
    end
    title('Aggregate Estimated Vibration PSD Overlay (Rigid Body assumption)');
    xlabel('Frequency [Hz]'); ylabel('Vibration [g^2/Hz]');
    legend(plot_handles_vib, legend_names_vib, 'Location', 'best');
    xlim([0, psd_data(1).f(end)]);
    hold off;
end