%% NASA CEA Relation Plotter
% Author: Andrew Radulovich (aradulov@purdue.edu), Kamon Blong
% First Created: 4/10/2023
% Last Updated: 11/30/2024

%{ 
Description: Plots ISP vs OF ratio of fuels for a given oxidizer
%}

clear; clc; close all;

%% Initialize Variables

% propellants
fuels = [{'CH3OH'}, {'C2H5OH'}, {'C3H8O,2propanol'}, {'Jet-A(L)'}];  % list of fuel strings for NASA CEA
%fuels = {'CH3OH'};  % list of fuel strings for NASA CEA
oxidizer = {'O2(L)'}; % oxidizer formula for NASA CEA
fuel_weight = 100;      % fuel weights
oxidizer_weight = 100;  % Oxidizer weight
fuel_temps = [273.15, 273.15, 273.15, 273.15];   % fuel temperatures [K]
oxidizer_temp = 90.19;   % oxidizer temperature [K]

% Pressures
P_c = 250;            % chamber pressure [psi]
P_e = 17;           % exit pressure [psi]
PcPe = P_c/P_e; % Ratio of chamber pressure to exit pressure

% Initialize OF and ISP matrix
min_value_OF = .5;
max_value_OF = 3;
step_value_OF = .04;

OF_matrix = min_value_OF:step_value_OF:max_value_OF;
isp_matrix = ones(length(fuels), round((max_value_OF - min_value_OF) / step_value_OF));

% miscellaneous
legend_str = {'Data 1', 'Data 2'};

%% CEA Calculations
j = 1;
% frozen calculations
for fuel = fuels
    CEA_out1 = callCEA('fr', P_c, 'psi', 'o/f',OF_matrix,'pip', PcPe, fuel, 'K', fuel_temps(j), fuel_weight, oxidizer, 'K', oxidizer_temp, oxidizer_weight  );
    isp_parse = squeeze(CEA_out1('isp'));
    isp_vector = isp_parse(:,2);
    isp_matrix(j,:) = transpose(isp_vector);
    propellant_string = num2cell(char(fuel + " / " + oxidizer + " frozen reactions"), 2);
    legend_str{j} = [strjoin(propellant_string), j];
    j = j + 1;
end
% equillibrium calculations
for fuel = fuels
    CEA_out2 = callCEA('eq', P_c, 'psi', 'o/f',OF_matrix,'pip', PcPe, fuel, 'K', fuel_temps(j-4), fuel_weight, oxidizer, 'K', oxidizer_temp, oxidizer_weight  );
    isp_parse = squeeze(CEA_out2('isp'));
    isp_vector = isp_parse(:,2);
    isp_matrix(j,:) = transpose(isp_vector);
    propellant_string = num2cell(char(fuel + " / " + oxidizer + ' equilibrium reactions'), 2);
    legend_str{j} = [strjoin(propellant_string), j];
    j = j + 1;
end

%% Plotting
figure('Name', 'Fuel Trade Study')
hold("on")
legend(legend_str, 'location', 'Northeast');
plot(OF_matrix, isp_matrix, 'Linewidth', 2.5);
%plot(OF_matrix, isp_matrix(1,:), OF_matrix, isp_matrix(2,:), OF_matrix, isp_matrix(3,:), OF_matrix, isp_matrix(4,:), OF_matrix, isp_matrix(5,:), OF_matrix, isp_matrix(6,:), 'Linewidth', 2.5);
legend(legend_str);
xlabel('Mixture Ratio')
ylabel('Ideal Isp (sec)')
%set(gca, 'XLim', [0.75, 2.4], 'FontSize', 17)
%set(gca, 'YLim', [205, 265])
grid on
title("Propellant Performance vs Mixture Ratio: " + P_c + " psi P_c, " + P_e + " psi P_e")