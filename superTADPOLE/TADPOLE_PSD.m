% Hotfire Detection and PSD Analysis Script
% Docks figures, fixes pre-ignition clipping, and overlays PSDs with confidence bounds.

clear; clc; close all;

% 1. Set figures to dock by default
set(0, 'DefaultFigureWindowStyle', 'docked');

% --- Configuration Parameters ---
dataFolder = 'Raw Data Logs';          % Set to 'Raw Data Logs' if script is outside the folder
filePattern = fullfile(dataFolder, 'TEST_*_data.mat');
testFiles = dir(filePattern);

if isempty(testFiles)
    error('No test data files found. Please check your directory path.');
end

% Detection & Clipping Thresholds
pc_min_operating = 100;   % [psi] Strict threshold: Must hit this to be a hotfire
envelope_thresh  = 50;    % [psi] Boundary threshold: Defines the start/end to survive throttles
pad_time         = 8;     % [sec] Time buffer for the VISUAL plot only

% Data storage for the final aggregate plot
psd_data = struct('name', {}, 'f', {}, 'pxx', {}, 'pxxc', {});

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
    fprintf('Hotfire Detected! Processing...\n');
    
    Fs = 1 / mean(diff(t));
    
    % 2. Improved Logic: Envelope Detection
    % Find the overall envelope of the burn using the lower boundary threshold.
    % This ensures purges (~20psi) are ignored, but deep throttles (~100psi)
    % won't accidentally break the continuity.
    envelope_indices = find(PC > envelope_thresh);
    if isempty(envelope_indices)
        continue; 
    end
    
    % Span from the very first crossing of the envelope to the very last
    best_start = envelope_indices(1);
    best_end = envelope_indices(end);
    
    % Extract strict steady-state region for PSD (trimming 0.25s off ends to skip startup/shutdown transients)
    trim_idx = round(0.25 * Fs);
    idx_psd_start = min(best_start + trim_idx, best_end);
    idx_psd_end   = max(best_end - trim_idx, best_start);
    
    PC_burn = PC(idx_psd_start:idx_psd_end);
    PC_zero_mean = PC_burn - mean(PC_burn);
    
    % 3. Calculate PSD with Confidence Bounds (95%)
    window = round(Fs * 0.25);
    noverlap = round(window * 0.5);
    [pxx, f, pxxc] = pwelch(PC_zero_mean, window, noverlap, [], Fs, 'ConfidenceLevel', 0.95);
    
    % Store for final plot
    psd_data(end+1).name = fileName;
    psd_data(end).f = f;
    psd_data(end).pxx = pxx;
    psd_data(end).pxxc = pxxc;
    
    % 4. Visual Plotting Extraction (Wider view for context)
    % We use burn_indices(1) to burn_indices(end) to pad around the whole event
    pad_idx = round(pad_time * Fs);
    idx_plot_start = max(1, envelope_indices(1) - pad_idx);
    idx_plot_end = min(length(t), envelope_indices(end) + pad_idx);
    
    t_clip  = t(idx_plot_start:idx_plot_end);
    PC_clip = PC(idx_plot_start:idx_plot_end);
    
    % 5. Individual Figure Generation
    figure('Name', sprintf('Analysis: %s', fileName));
    
    % Subplot A: Clipped Autosequence Time History
    subplot(1, 2, 1);
    plot(t_clip, PC_clip, 'r', 'LineWidth', 1.2);
    grid on; hold on;
    title(sprintf('%s - Clipped Autosequence', strrep(fileName, '_', '\_')));
    xlabel('Time [sec]'); ylabel('Chamber Pressure [psi]');
    
    % Highlight the exact region used for the PSD
    xline(t(idx_psd_start), '--k', 'PSD Start', 'LabelOrientation', 'horizontal');
    xline(t(idx_psd_end), '--k', 'PSD End', 'LabelOrientation', 'horizontal');
    xlim([min(t_clip) max(t_clip)]);
    
    % Subplot B: Frequency Domain (PSD)
    subplot(1, 2, 2);
    % Plot confidence bounds (convert to dB)
    f_fill = [f; flipud(f)];
    p_fill = 10*log10([pxxc(:,1); flipud(pxxc(:,2))]);
    fill(f_fill, p_fill, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
    % Plot main PSD line
    plot(f, 10*log10(pxx), 'b', 'LineWidth', 1.2); 
    grid on;
    title('PSD - Steady State Burn Region');
    xlabel('Frequency [Hz]'); ylabel('Power/Frequency [dB/Hz]');
    xlim([0, Fs/2]);
end

fprintf('\nBatch processing complete. Generating overlay plot...\n');

% 6. Final Aggregate Overlay Figure
if ~isempty(psd_data)
    figure('Name', 'Aggregate PSD Overlay');
    hold on; grid on;
    
    color_map = lines(length(psd_data)); % Generate distinct colors
    plot_handles = gobjects(length(psd_data), 1);
    legend_names = cell(length(psd_data), 1);
    
    for k = 1:length(psd_data)
        f_k = psd_data(k).f;
        pxx_k = 10*log10(psd_data(k).pxx);
        pxxc_k = 10*log10(psd_data(k).pxxc);
        
        % Plot transparent confidence bounds
        f_fill = [f_k; flipud(f_k)];
        p_fill = [pxxc_k(:,1); flipud(pxxc_k(:,2))];
        fill(f_fill, p_fill, color_map(k,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
        % Plot main line and store handle for legend
        plot_handles(k) = plot(f_k, pxx_k, 'Color', color_map(k,:), 'LineWidth', 1.5);
        legend_names{k} = strrep(psd_data(k).name, '_', '\_');
    end
    
    title('Aggregate PSD Overlay (with 95% Confidence Bounds)');
    xlabel('Frequency [Hz]');
    ylabel('Power/Frequency [dB/Hz]');
    legend(plot_handles, legend_names, 'Location', 'best');
    xlim([0, psd_data(1).f(end)]);
    hold off;
end