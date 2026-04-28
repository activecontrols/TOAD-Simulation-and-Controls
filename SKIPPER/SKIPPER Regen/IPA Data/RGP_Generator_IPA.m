%% RGP_Full_Generator_IPA.m
% 1. SETUP AXES
T_range = linspace(280, 350, 20); % Extended to include boiling point (~355K)
P_range = linspace(1e5, 20e6, 25); % Pa
[TT, PP] = meshgrid(T_range, P_range);
n_T = length(T_range); n_P = length(P_range);

% 2. ALLOCATE ALL TABLES
% Liquid (101-106) and Vapor (1-6)
L_tables = cell(1,6); V_tables = cell(1,6);
for i=1:6, L_tables{i}=zeros(n_P,n_T); V_tables{i}=zeros(n_P,n_T); end

% 3. POPULATE DATA
R_ipa = 8314 / 60.10; % Gas constant for IPA (J/kgK)

for i = 1:n_P
    for j = 1:n_T
        % --- LIQUID DATA (Using your function) ---
        [rho, cp, visc, cond] = IPA_properties(T_range(j), P_range(i)/1000);
        L_tables{3}(i,j) = 1/rho;              % Spec Volume (103)
        L_tables{4}(i,j) = cp * 1000;          % Cp J/kg (104)
        L_tables{1}(i,j) = L_tables{4}(i,j) * (T_range(j) - 273.15); % Enthalpy (101)
        L_tables{2}(i,j) = 1140;               % Speed of Sound (102)
        L_tables{5}(i,j) = visc;               % Viscosity (105)
        L_tables{6}(i,j) = cond;               % Conductivity (106)
        
        % --- VAPOR DATA (Ideal Gas Placeholders for File Integrity) ---
        V_tables{3}(i,j) = (R_ipa * T_range(j)) / P_range(i); % v = RT/P (3)
        V_tables{4}(i,j) = 1500;               % Generic Cp for vapor (4)
        V_tables{1}(i,j) = V_tables{4}(i,j) * (T_range(j) - 273.15); % h (1)
        V_tables{2}(i,j) = sqrt(1.2 * R_ipa * T_range(j)); % a (2)
        V_tables{5}(i,j) = 1e-5;               % mu (5)
        V_tables{6}(i,j) = 0.02;               % k (6)
    end
end

% 4. WRITE THE COMPLETE FILE
fid = fopen('IPA_Complete.rgp', 'w');
fprintf(fid, '$$$$DATA\n$$$IPA_LIQUID\n  NAME IPA_Liquid\n  MOLAR_MASS 60.10\n');
fprintf(fid, '  CRITICAL_P 4.763e6\n  CRITICAL_T 508.3\n');

% Write Vapor Tables (1-6)
write_rgp_block(fid, V_tables, [1:6], P_range, T_range);
% Write Liquid Tables (101-106)
write_rgp_block(fid, L_tables, [101:106], P_range, T_range);

fclose(fid);

function write_rgp_block(fid, data_cell, ids, P, T)
    for k = 1:6
        fprintf(fid, '$TABLE %d\n', ids(k));
        fprintf(fid, '  P_MIN %e\n  P_MAX %e\n', min(P), max(P));
        fprintf(fid, '  T_MIN %e\n  T_MAX %e\n', min(T), max(T));
        fprintf(fid, '  N_P %d\n  N_T %d\n', length(P), length(T));
        for r = 1:length(P)
            fprintf(fid, '  %e', data_cell{k}(r, :));
            fprintf(fid, '\n');
        end
    end
end