% ********************************************************************
%
% PSP Throttle Control
% File Name: Post_sim_data_tables.m
% Author: Eric Umminger
% Date: 9/26/26
%
% Program Description: This script extracts the workspace data from the
% Simulink simulation and creates the data tables.
%
% ********************************************************************

% ********************************************************************
% Extracting timeseries data

cv_ox_cmd_timeseries = out.cv_ox_cmd;
pc_cmd_timeseries = out.pc_cmd;

cv_fu_cmd_timeseries = out.cv_fu_cmd;
pox_i_m_timeseries = out.pox_i_m;
pfu_i_cmd_timeseries = out.pfu_i_cmd;

% ********************************************************************
% Extracting matrices

cv_ox_cmd = cv_ox_cmd_timeseries.Data(:);
pc_cmd = pc_cmd_timeseries.Data(:);

cv_fu_cmd = cv_fu_cmd_timeseries.Data(:);
pox_i_m = pox_i_m_timeseries.Data(:);
pfu_i_cmd = pfu_i_cmd_timeseries.Data(:);

% ********************************************************************
% Creating tables

pcCmd_to_cvOxCmd = sortrows([pc_cmd, cv_ox_cmd], 1, 'ascend');
poxIM_to_pfuICmd = sortrows([pox_i_m, pfu_i_cmd], 1, 'ascend');
pfuICmd_to_cvFuCmd = sortrows([pfu_i_cmd, cv_fu_cmd], 1, 'ascend');

% ********************************************************************
% Writing to files

% Chamber pressure to oxygen valve coefficient
T1 = array2table(pcCmd_to_cvOxCmd, 'VariableNames', {'pc_cmd', 'cv_ox_cmd'});
writetable(T1, 'pcCmd_to_cvOxCmd.txt', 'Delimiter', '\t');

% Measured oxygen injector pressure to commanded fuel injector pressure
T2 = array2table(poxIM_to_pfuICmd, 'VariableNames', {'pox_i_m', 'pfu_i_cmd'});
writetable(T2, 'poxIM_to_pfuICmd.txt', 'Delimiter', '\t');

% Commanded fuel injector pressure to commanded fuel valve coefficient
T3 = array2table(pfuICmd_to_cvFuCmd, 'VariableNames', {'pfu_i_cmd', 'cv_fu_cmd'});
writetable(T3, 'pfuICmd_to_cvFuCmd.txt', 'Delimiter', '\t');

fprintf("\nTables created.\n\n")