%% NASA CEA Relation Plotter
% Authors: Andrew Radulovich (aradulov@purdue.edu), Kamon Blong
% First Created: 4/10/2023
% Last Updated: 11/30/2024

%{ 
Description: Plots ISP vs OF ratio of oxidizers for a given fuel
%}

clear; clc; close all;

%% Initialize Variables

% propellants
fuel = {'C3H8O,2propanol'};  % list of fuel strings for NASA CEA
%fuel = {'CH3OH'};
oxidizers = [{'O2(L)'}, {'H2O2(L)'}, {'N2O'}]; % oxidizer formula for NASA CEA
fuel_weight = 100;      % fuel weights
oxidizer_weight = 100;  % Oxidizer weight
fuel_temp = 273.15;   % fuel temperature [K]
oxidizer_temp = [90.19, 273.15, 273.15];   % fuel temperature [K]

% Pressures
P_c = 250;            % chamber pressure [psi]
P_e = 17;           % exit pressure [psi]
PcPe = P_c/P_e; % Ratio of chamber pressure to exit pressure

% Initialize OF and ISP Matrix
min_value_OF = .5;
max_value_OF = 7;
step_value_OF = .1;

OF_matrix = min_value_OF:step_value_OF:max_value_OF;
isp_matrix = ones(length(oxidizers), 1 + round((max_value_OF - min_value_OF) / step_value_OF));

% miscellaneous
legend_str = {'Data 1', 'Data 2'};

%% CEA Calculations
j = 1;

for oxidizer = oxidizers
    CEA_out0 = callCEA('fr', P_c, 'psi', 'o/f',OF_matrix,'pip', PcPe, fuel, 'K', fuel_temp, fuel_weight, oxidizer, 'K', oxidizer_temp(j), oxidizer_weight  );
    isp_parse = squeeze(CEA_out0('isp'));
    isp_vector = isp_parse(:,2);
    isp_matrix(j,:) = transpose(isp_vector);
    propellant_string = num2cell(char(fuel + " / " + oxidizer), 2);
    legend_str{j} = [strjoin(propellant_string), j];
    j = j + 1;
end

%% Plotting
figure('Name', 'Fuel Trade Study')
hold("on")
legend(legend_str, 'location', 'Northeast');
plot(OF_matrix, isp_matrix, 'Linewidth', 2.5);
legend(legend_str);
xlabel('Mixture Ratio')
ylabel('Ideal Isp (sec)')
set(gca, 'XLim', [0.75, 5.25], 'FontSize', 17)
set(gca, 'YLim', [175, 255])
grid on
title("Propellant Performance vs Mixture Ratio: " + P_c + " psi P_c, " + P_e + " psi P_e")