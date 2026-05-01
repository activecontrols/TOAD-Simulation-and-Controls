% Regenerative Cooling Sizing Code
% Authors: Adam Grendys
% Start Date: 2/27/26
% Description: This code is based off of PSP:AC's "size_regen" script written by Grant Williams, Zach Hodgdon, Andrew Radulovich, Alex Suppiah, Jan Ayala, Kamon Blong. 

function [Lifespan, PressDrop] = SKIPPERRegen(Data, NumChannels, WallThickness, ChannelHeight, ChannelWidth, DisplayMode)
New_CEA = false;
fclose all;
close all;
u = convertUnits;
CEA_input_name = 'regrendysCEA';
dfTol = 5e-6;
optionsMoody = optimset('TolX', dfTol, 'Display', 'off');
tic

%% SIMULATION PARAMETERS
steps = 100; % Number of steps along chamber (Change resolution of simulation)
materialchoice = 4; % 0: 6061-RAM2, 1: Inconel718, 2: GrCop42, 3: Inconel625, 4: C101/OFHC Copper, 5: AlSi10Mg
heat_correlation = 1; % 1 = Bartz, 2 = Heister Bartz
coolant_direction = 0; % 0 = opposite flow, 1 = flow with hot gas 
FEA_outputs = 0; % 1 = yes, 0 = no
dogleg = 0; % 1 = yes, 0 = no, supertadpole regen channel dogleg at injector
traditional = 1;  % 1 = yes, 0 = no, changes how channel dimensions are interpolated for a traditonal vs. printed chamber

throttle = 0.5; % throttle percent - e.g. 1 = 100%, 0.5 = 50%
num_channels = round(NumChannels); % number of regenerative cooling channels      
coolant = "isopropyl alcohol"; % coolant definition ("isopropyl alcohol", "water", "methanol", "ethanol")
fuel = {'C3H8O,2propanol'}; % fuel definition
oxidizer = 'O2(L)'; % oxidizer definition

T_amb = 297; % Atmospheric temp [K]                            
fuel_temp = throttle * -25.556 + 427.605; % [K], linear fit for predicted regen outlet temps
oxidizer_temp = 90.17; % [K]

P_c = throttle * 250; % chamber pressure [psi] 
P_e = throttle * 16.5; % exit pressure [psi]
P_inlet = 375 * throttle + 200; % Regen inlet pressure [psi]  
total_OF = 1.2; % oxidizer/fuel ratio (TOTAL, including film)  
total_mdot = throttle * 3.4178026327 / 2.205; % TOTAL chamber mass flow [kg/s]  
mdot_coolant = total_mdot / (1 + total_OF); % Coolant/fuel mass flow [kg/s]

qdot_tolerance = 0.0001; % Heat loop convergence criteria
debug = 0; % Debug tool (1 = on, 0 = off)

% Channel Defintion
t_w = WallThickness .* 0.0254; % channel wall thickness [1 min 2] [m] 
h_c = ChannelHeight .* 0.0254; % channel height [1 min 2] [m]   
w_c = ChannelWidth .* 0.0254; % channel width [1 min] [m]  

heatflux_factor = -0.10833 * throttle + 0.6433; % Scaling factor [0 to 1], Linear Fit to Tadpole Data 

%% CALCULATIONS
%----------------------------------%

%% Material Properties
% if materialchoice == 0
%     properties = readmatrix(pwd + "/Material Data/al6061_RAM2.xlsx");
%     k_w = properties(16:24,1:2); % thermal conductivity [W/m-K]
%     E = [properties(16:21, 8) properties(16:21,9)];
%     CTE = [properties(16:20,4) properties(16:20,5)]; % [m/m*K]
%     yield_strength = properties(1:11,1:2);
%     elongation_break = [properties(1:11,1) properties(1:11,4)];
% elseif materialchoice == 1 
%     properties = readmatrix(pwd + "/Material Data/Inconel718.xlsx");
%     k_w = properties(13:end,1:2); % thermal conductivity [W/m-K]
%     yield_strength = properties(1:7,1:2);
%     elongation_break = [properties(1:7,1) properties(1:7,5)];
%     E = [properties(1:6, 9) properties(1:6,10)]; %youngs modulus
%     CTE = [properties(1:7,1) properties(1:7,3)]; % [ppm]
% elseif materialchoice == 2
%     properties = readmatrix(pwd + "/Material Data/GrCop42.xlsx");
%     k_w = properties(13:end,1:2); % thermal conductivity [W/m-K] 
%     E = [properties(1:13, 9) properties(1:13,10)];
%     CTE = [properties(1:5,1) properties(1:5,3)]; % [ppm] 
%     yield_strength = properties(1:8,1:2);
%     elongation_break = [properties(1:8,1) properties(1:8,5)];
% elseif materialchoice == 3
%     properties = readmatrix(pwd + "/Material Data/Inconel625.xlsx");
%     k_w = properties(8:19,1:2); % thermal conductivity [W/m-K] 
%     E = [properties(1:5, 1) properties(1:5,5)];
%     CTE = [properties(22:end,1:2)]; % [ppm] 
%     yield_strength = properties(1:5,1:2);
%     elongation_break = [properties(1:5,1) properties(1:5,3)];
% elseif materialchoice == 4
%     properties = readmatrix(pwd + "/Material Data/C101.xlsx");
%     k_w = properties(13:21,1:2); % thermal conductivity [W/m-K] 
%     E = [properties(1:5, 9:10)];
%     CTE = [properties(1:7,1) properties(1:7,3)]; % [ppm] 
%     yield_strength = properties(1:7, 1:2);
%     elongation_break = [properties(1:7,1) properties(1:7,5)];
% elseif materialchoice == 5
%     properties = readmatrix(pwd + "/Material Data/AlSi10Mg.xlsx");
%     k_w = properties(13:15,1:2); % thermal conductivity [W/m-K] 
%     E = [properties(1:6, 9:10)];
%     CTE = [properties(1:4,1) properties(1:4,3)]; % [ppm] 
%     yield_strength = properties(1:4,1:2);
%     elongation_break = [properties(1:4,1) properties(1:4,5)];
% end
%% Only props for copper
properties = Data.Properties;
k_w = properties(13:21,1:2); % thermal conductivity [W/m-K] 
E = [properties(1:5, 9:10)];
CTE = [properties(1:7,1) properties(1:7,3)]; % [ppm] 
yield_strength = properties(1:7, 1:2);
elongation_break = [properties(1:7,1) properties(1:7,5)];

v = 0.33; % GUESS poissons ratio
N = 35*8; % engine lifespan (for margin math)
SF = 4; % lifespan saftey factor (4 per NASA 5012C)

% Contour Interpolation
% contour = readmatrix('contour_SKIPPER_250_40.xlsx'); % import engine contour
contour = Data.Contour;
r_contour = contour(:,2)'; % contour radius [in]
x_contour = contour(:,1)'; % contour x-axis [in]
r_t = min(r_contour); % throat radius, throat location (in) 
total_length = abs(min(contour(:,1))) + max(contour(:,1)); % total length (in) 
chamber_length = contour(1,3);
converging_length = contour(2,3);
diverging_length = contour(3,3);

x_interpolated = linspace(min(x_contour), max(x_contour), steps);
r_interpolated = interp1(x_contour, r_contour, x_interpolated, 'linear', 'extrap'); % linearly interpolate radius vector  
dx = x_interpolated(2) - x_interpolated(1); % Small change in step length [in.]

% Channel Interpolation
t_w_x = NaN(1, steps);
w_c_x = NaN(1, steps);
h_c_x = NaN(1, steps);

if traditional 
    % Channel Interior Surface Roughness
    roughness_abs = (125 * 0.0254) * 10^-6; % Surface roughness [m]
    
    PressureLoss_factor = 1.5; % Scaling factor for coolant pressure loss (to be verified)
    P_minor_coefs = [2, 0.5, 1];  % Minor Loss Coefficients [Inlet Cv Manifold, Throat Bend, Injector Turnaround] 

    for i = 1:length(x_interpolated)
        if x_interpolated(i) <= -converging_length
            t_w_x(i) = t_w(1);
            w_c_x(i) = w_c(1);
            h_c_x(i) = h_c(1);
        elseif x_interpolated(i) > -converging_length && x_interpolated(i) <= 0
            t_w_x(i) = (t_w(2) - t_w(1)) / converging_length * x_interpolated(i) + t_w(2);
            w_c_x(i) = w_c(2);
            h_c_x(i) = h_c(2);
        else
            t_w_x(i) = (t_w(3) - t_w(2)) / diverging_length * x_interpolated(i) + t_w(2);
            w_c_x(i) = w_c(2);
            h_c_x(i) = (h_c(3) - h_c(2)) / diverging_length * x_interpolated(i) + h_c(2);
        end
    end
else
    % Channel Interior Surface Roughness
    roughness_table = readmatrix(pwd + "/Material Data/surface_roughness.xlsx",'Range','A20:B24'); % High-End Roughness
    roughness_abs = roughness_table(2,2) * 10^-6; % Surface roughness [m] [45, 90]
    
    
    PressureLoss_factor = 5; % Scaling factor for coolant pressure loss (fit to tadpole data)
    P_minor_coefs = [4.5, 1, 2];  % Minor Loss Coefficients [Inlet Cv Manifold, Throat Bend, Injector Turnaround]

    for i = 1:length(x_interpolated)
        if x_interpolated(i) <= -converging_length
            t_w_x(i) = t_w(1);
            w_c_x(i) = w_c(1);
            h_c_x(i) = h_c(1);
        elseif x_interpolated(i) > -converging_length && x_interpolated(i) <= 0
            t_w_x(i) = (t_w(2) - t_w(1)) / converging_length * x_interpolated(i) + t_w(2);
            w_c_x(i) = (w_c(2) - w_c(1)) / converging_length * x_interpolated(i) + w_c(2);
            h_c_x(i) = (h_c(2) - h_c(1)) / converging_length * x_interpolated(i) + h_c(2);
        else
            t_w_x(i) = (t_w(3) - t_w(2)) / diverging_length * x_interpolated(i) + t_w(2);
            w_c_x(i) = (w_c(3) - w_c(2)) / diverging_length * x_interpolated(i) + w_c(2);
            h_c_x(i) = (h_c(3) - h_c(2)) / diverging_length * x_interpolated(i) + h_c(2);
        end
    end
end

% Unit Control
P_c = P_c * 6894.76; % Chamber presssure (psi to Pa)
P_e = P_e * 6894.76; % Exit Pressure (psi to Pa)
P_inlet = P_inlet * 6894.76; % Regen inlet pressure (psi to Pa)
x_interpolated = x_interpolated * 0.0254; % x-vector (in to m)
r_interpolated = r_interpolated * 0.0254; % radius-vector (in to m)
r_t = r_t * 0.0254; % throat radius (in to m)
total_length = total_length * 0.0254; % length (in to m)
chamber_length = chamber_length * 0.0254; % length (in to m)
converging_length = converging_length * 0.0254; % length (in to m)
diverging_length = diverging_length * 0.0254; % length (in to m)
dx = dx * 0.0254; % step length (in to m)

if dogleg
    l_straight = 0.442 * 0.0254; % lenght of straight section with highest wall thickness (m)
    max_wt = 0.4768 * 0.0254; % wall thickness at farthest point (m)
    leg_angle = 45; % dogleg angle (deg)

    l_leg = (max_wt - t_w(1)) / sind(leg_angle); % axial length of dogleg (m)
    index_straight = floor(l_straight / dx);
    index_leg = floor(l_leg / dx);

    t_w_x(1:index_straight) = max_wt;
    
    for j = 1:index_leg
        t_w_x(j + index_straight) = max_wt - j * (max_wt - t_w(1)) / index_leg;
        if t_w_x(j + index_straight) < 0
            t_w_x(j + index_straight) = t_w(1);
        end
    end
end

% Channel centerline arc length per axial station
r_channel_center = r_interpolated + t_w_x + 0.5 * h_c_x;

ds_vec = zeros(1, steps);
for k = 1:steps-1
    ds_vec(k) = sqrt((x_interpolated(k+1) - x_interpolated(k))^2 + (r_channel_center(k+1) - r_channel_center(k))^2 );
end
ds_vec(end) = ds_vec(end-1);

%% CHAMBER HEAT TRANSFER
% Step 1: Initialize Properties 
mdot_channel = mdot_coolant / num_channels; % per channel mdot [kg/s]
A_t = pi * r_t^2; % throat area (m^2)
D_t = 2 * r_t; % Diameter of throat (m)
R_of_curve = 1.5 * r_t; % Radius of Curvature at throat (input)

Qtot = 0; %total heat through entire chamber [W]

% Area Ratios
subsonic_area_ratios = (pi * r_interpolated(x_interpolated < 0) .^ 2) / A_t; % subsonic area ratios
supersonic_area_ratios = (pi * r_interpolated(x_interpolated >= 0) .^ 2) / A_t; %  supersonic area ratios (including throat)
A_ratio = (pi * r_interpolated.^2) / A_t;

% Minor Pressure Loss Locations
idx_Cv = steps - 1; % Index of Cv inlet manifold
[~, idx_throat] = min(r_interpolated); % Index of throat
idx_inj = 2; % Index of injector turnaround

% Fluid Property Vectors
P_coolant = zeros(1, steps); % coolant pressure
T_r = zeros(1, steps); % recovery temperature
T_coolant = zeros(1, steps); % coolant temp
rho_coolant = zeros(1, steps); % coolant density
vel_coolant = zeros(1, steps); % coolant velocity
cp_coolant = zeros(1, steps); % coolant specific heat (constant presure)
mu_coolant = zeros(1, steps); % coolant dynamic visocosity
k_coolant = zeros(1, steps); % coolant thermal conductivity
vel_g = zeros(1, steps); % chamber velocity
Re_coolant = zeros(1, steps); % coolant Reynolds number
Pr_coolant = zeros(1, steps); % coolant Prandtl number
boilingPt = zeros(1,steps); %coolant boiling point
Qdot_l = zeros(1, steps);  % liquid convective heat transfer
Qdot_g = zeros(1, steps);  % gas convective heat transfer
heatflux = zeros(1, steps); % Heat Flux [W/m^2]
heatflux_fin = zeros(1, steps); % Fin Heat Flux [W/m^2]
T_wl = zeros(1, steps);    % liquid wall temperature
T_wg = zeros(1, steps);    % gas wall temperature
h_g = zeros(1, steps);     % gas film coefficient
sigma = zeros(1, steps);   % film coefficient correction factor
h_l = zeros(1, steps);     % liquid film coefficient
Nu_l = zeros(1,steps);     % Nusselt Number

% Material Property Vectors
k_w_current = zeros(1,steps);
yield = zeros(1,steps);
elong = zeros(1,steps);
E_current = zeros(1,steps);
CTE_current = zeros(1,steps);
CTE_liq_side = zeros(1,steps);

% Fin Vectors
fin_w = zeros(1,steps);
A_fin = zeros(1,steps);
Qdot_fin = zeros(1,steps);
eta_fin = zeros(1,steps);

% Stress Vectors
epsilon_emax = zeros(1,steps);
sigma_t = zeros(1,steps); % tangential stress
sigma_tp = zeros(1,steps); % tangential stress pressure
sigma_tt = zeros(1,steps); % tangential stress temp
sigma_l = zeros(1,steps); % longitudinal stress
sigma_ll = zeros(1,steps); % longitudinal stress
sigma_lc = zeros(1,steps); % longitudinal stress
sigmab = zeros(1,steps); % buckling stress
sigma_v = zeros(1,steps); % von mises stress
sigma_vl = zeros(1,steps); % von mises stress
sigma_vc = zeros(1,steps); % von mises stress
sigma_tp_cold = zeros(1,steps); % Pressing channels before hotfire
epsilon_lc = zeros(1,steps);
epsilon_ll = zeros(1,steps); 
epsilon_t = zeros(1,steps);
epsilon_vc = zeros(1,steps);
epsilon_vl = zeros(1,steps);
epsilon_tp = zeros(1,steps);
epsilon_tt = zeros(1,steps);

epsilon_tota = zeros(1,steps);
epsilon_tott = zeros(1,steps);
epsilon_pa = zeros(1,steps);
epsilon_pt = zeros(1,steps);
epsilon_peff = zeros(1,steps);
epsilon_cs = zeros(1,steps);
MS = zeros(1,steps);
num_fires = zeros(1,steps);
num_fires_lowcycle = zeros(1,steps);
epsilon_toteff = zeros(1,steps);
epsilon_emaxeff = zeros(1,steps);
epsilon_eeff = zeros(1,steps);
sigma_eff = zeros(1,steps);
sigma_a = zeros(1,steps);
sigma_t2 = zeros(1,steps);
epsilon_cs_tot = zeros(1,steps);
epsilon_cs_spacex = zeros(1,steps);
MS_spacex = zeros(1,steps);
MS_lowcycle = zeros(1,steps);
deltaT1 = zeros(1,steps);
deltaT2 = zeros(1,steps);

if debug
    Qdot_check = zeros(1, steps);
end

counter = 1;
error = zeros(1, steps); % Heat loop convergence error

%% Call CEA for all Area Ratios
if New_CEA
    fprintf("Running CEA. %d Iterations to run. (~%.2f minutes) \nIteration: ",steps,((steps*1.48-14)/60)); % gives time estimate to running CEA
    i = 1;
    for sub = subsonic_area_ratios
        [~, ~, ~, M(i), gamma(i), P_g(i), T_g(i), rho_g(i), mu_g(i), Pr_g(i), ~, ~, ~, cp_g(i)] = RunCEA(P_c, P_e, fuel, 100, fuel_temp, oxidizer, oxidizer_temp, total_OF, sub, 0, 2, 0, 0, CEA_input_name);
        fprintf("%d ",i);
        i = i + 1;
    end
    i = size(subsonic_area_ratios, 2) + 1;
    for sup = supersonic_area_ratios
        [~, ~, ~, M(i), gamma(i), P_g(i), T_g(i), rho_g(i), mu_g(i), Pr_g(i), ~, ~, ~, cp_g(i)] = RunCEA(P_c, P_e, fuel, 100, fuel_temp, oxidizer, oxidizer_temp, total_OF, 0, sup, 2, 0, 0, CEA_input_name);
        fprintf("%d ",i);
        i = i + 1;
    end
    
    [c_star, ~, ~, ~, ~, P_g_tot, T_g_tot, rho_g_tot, mu_g_tot, Pr_g_tot, ~, ~, ~, cp_g_tot] = RunCEA(P_c, P_e, fuel, 100, fuel_temp, oxidizer, oxidizer_temp, total_OF, 0, 0, 1, 0, 0, CEA_input_name);
    
    % Save CEA data to .mat file
    save("RegenCEA_Data.mat","M","gamma","P_g","T_g","rho_g","mu_g","Pr_g","cp_g","c_star","P_g_tot","T_g_tot","rho_g_tot","mu_g_tot","Pr_g_tot","cp_g_tot"); 
else
    % Load previously calculated CEA data
    RegenCEA_Data = load("RegenCEA_Data.mat");
    M = RegenCEA_Data.M; gamma = RegenCEA_Data.gamma;
    P_g = RegenCEA_Data.P_g; T_g = RegenCEA_Data.T_g;
    rho_g = RegenCEA_Data.rho_g; mu_g = RegenCEA_Data.mu_g;
    Pr_g = RegenCEA_Data.Pr_g; cp_g = RegenCEA_Data.cp_g;
    c_star = RegenCEA_Data.c_star; P_g_tot = RegenCEA_Data.P_g_tot;
    T_g_tot = RegenCEA_Data.T_g_tot; rho_g_tot = RegenCEA_Data.rho_g_tot;
    mu_g_tot = RegenCEA_Data.mu_g_tot; Pr_g_tot = RegenCEA_Data.Pr_g_tot;
    cp_g_tot = RegenCEA_Data.cp_g_tot;
end

% Boundary Conditions
if coolant_direction 
    P_coolant(1) = P_inlet;
    T_coolant(1) = T_amb;

    points = 1:steps;
elseif ~coolant_direction 
    P_coolant(steps) = P_inlet;
    T_coolant(steps) = T_amb;

    points = 1:steps;
    points = steps + 1 - points; %reverses direction - i.e from [1 2 3] to [3 2 1]
else
    fprintf("/nCHOOSE VALID COOLANT DIRECTION")
end

% Initial Temperature Guesses
r = Pr_g.^ (1 / 3); % Recovery factor
T_wg = T_g .* (1 + (gamma - 1) / 2 .* r .* M .^ 2); % hotwall guess = recovery temp [K]
T_wl = T_coolant; % coldwall guess = coolant start temp [K]

for i = points % 1 = injector, steps = exit
    T_wg_min = T_coolant(i); % minimum temperature bound [K]
    T_wg_max = 3000; % maximum temperature bound [K]

    % Fluid Area
    A_channel = h_c_x(i) * w_c_x(i);
    p_wet = 2 * h_c_x(i) + 2 * w_c_x(i);
    hydraulic_D = 4 * A_channel / p_wet; % channel hydrualic diameter [m]

    counter = 1; % Heat Convergence Loop counter
    converged = 0; % Heat Convergence Loop end parameter

    while ~converged
        deltaP_minor = 0;
        r = Pr_g(i) ^ (1 / 3); % recovery factor for a turbulent free boundary layer [N/A] - biased towards larger engines, very small engines should use Pr^.5 (Heister Table 6.2).
        T_r(i) = T_g(i) * (1 + (gamma(i) - 1) / 2 * r * M(i) ^ 2); % recovery temperature [K] - corrects for compressible boundry layers (Heister EQ 6.15). 

        % Heat Transfer Correlation Selection
        if heat_correlation == 1
            % Bartz Correlation
            w = log(mu_g(i) / mu_g_tot) / log(T_g(i) / T_g_tot); if isnan(w); w = 0; end
            sigma(i) = ((0.5 * (T_wg(i) / T_g_tot) * (1 + ((gamma(i) - 1) / 2) * M(i)^2) + 0.5) ^ (0.8 - w/5) * (1 + (gamma(i) - 1) / 2 * M(i)^2) ^ (w/5))^(-1); % film coefficient correction factor [N/A] (Huzel & Huang eq. 4-14 (pg.86))
            h_g(i) = heatflux_factor*(0.026 / (D_t ^ 0.2)) * ((mu_g(i) ^ 0.2) * cp_g(i) / (Pr_g(i) ^ 0.6)) * ((P_g_tot / c_star) ^ 0.8) * ((D_t / R_of_curve) ^ 0.1) * ((1/A_ratio(i)) ^ 0.9) * sigma(i); % gas film coefficient [W/m^2-K] - bartz equation (Huzel & Huang pg.86).

        elseif heat_correlation == 2
            % Revised Bartz Correlation
            w = log(mu_g(i) / mu_g_tot) / log(T_g(i) / T_g_tot); if isnan(w); w = 0; end
            T_am = (T_wg(i) + T_g(i))/2;
            vel_g(i) = M(i) * sqrt(cp_g(i) * (gamma(i)-1) * T_r(i)); % velocity of the combustion flow, m/s
            h_g(i) = heatflux_factor * (0.026 / (2 * r_interpolated(i))^0.2) * (mu_g(i)^0.2 * cp_g(i) / (Pr_g(i) ^0.6)) * (rho_g(i) * vel_g(i))^0.8 * (T_g(i) / T_am)^0.8 * (T_am / T_g_tot)^(0.2 * w); % gas film coefficient [W/m^2-K] - revised bartz equation
        else
            fprintf("\n MUST SELECT VALID HEAT TRANSFER CORRELATION")
            converged = 1; % Loop End
        end

        % RESISTOR 1: Combustion Gas - Hotwall Convection
        A_hot = 2*pi * r_interpolated(i) / num_channels * ds_vec(i); % hotwall area [m^2]
        Qdot_g(i) = h_g(i) * A_hot * (T_r(i) - T_wg(i)); % gas convective heat transfer [W]

        % RESISTOR 2: Hotwall - Cold wall Conduction
        k_w_current(i) = interp1(k_w(:,1), k_w(:,2), (T_wg(i)+T_wl(i))/2, 'linear', 'extrap');
        T_wl(i) = T_wg(i) - (Qdot_g(i) / A_hot) * t_w_x(i) / k_w_current(i); % liquid wall temperature calculated via conduction through wall [K] (Heister EQ 6.29).        
        
        % RESISTOR 3: Coldwall - Coolant Convection
        if coolant == "isopropyl alcohol"
            P_kPa = P_coolant(i) / 1000;
            T_K = T_coolant(i);
        
            boilingPt(i)   = Data.F_boil(P_kPa);
            rho_coolant(i) = Data.F_rho(T_K, P_kPa);
            cp_coolant(i)  = Data.F_cp(T_K, P_kPa) * 1000; 
            mu_data = Data.F_mu_data(T_K);
            k_data  = Data.F_k_data(T_K);
        
            if (T_K - 273.15) < 36.5
                mu_fit = 4.5054 * exp(-0.031 * (T_K - 273.15)) / 1000;
            else
                mu_fit = 3.3724 * exp(-0.023 * (T_K - 273.15)) / 1000;
            end
            
            mu_avg = (mu_fit + mu_data) / 2;
            mu_coolant(i) = mu_avg * exp(5.4574e-9 * (P_coolant(i) - 101325));
            k_coolant(i)  = k_data * (1 + 2.9531e-9 * (P_coolant(i) - 101325));
        else
            mu_coolant(i) = double(py.CoolProp.CoolProp.PropsSI('V','T', T_coolant(i), 'P', P_coolant(i), coolant)); % viscosity of bulk coolant [Pa-s]
            cp_coolant(i) = double(py.CoolProp.CoolProp.PropsSI('C' , 'T', T_coolant(i), 'P', P_coolant(i), coolant)); % specific heat of coolant [J/kg-k] 
            k_coolant(i) = double(py.CoolProp.CoolProp.PropsSI('L', 'T', T_coolant(i), 'P', P_coolant(i), coolant)); % thermal conductivity of coolant [W/m-K]
            rho_coolant(i) = double(py.CoolProp.CoolProp.PropsSI('D','T', T_coolant(i),'P', P_coolant(i), coolant)); % density of the coolant [kg/m^3]
            boilingPt(i) = double(py.CoolProp.CoolProp.PropsSI('T','Q',0,'P', P_coolant(i), coolant)); % boiling point of the coolant 
        end

        vel_coolant(i) = mdot_channel / rho_coolant(i) / A_channel;
        Re_coolant(i) = rho_coolant(i) * vel_coolant(i) * hydraulic_D / mu_coolant(i); % reynolds number for channel flow [N/A] (Huzel and Huang , pg 90)
        Pr_coolant(i) = cp_coolant(i) * mu_coolant(i) / k_coolant(i); % prandtl number [N/A] (Huzel and Huang, pg 90) 

        f = moody(roughness_abs(1) / hydraulic_D, Re_coolant(i), optionsMoody); % friction factor
        % Accurate to a percent relative to fun moody.m call. Produces
        % ~10psi drop error. 
        % f = haaland(roughness_abs(1) / hydraulic_D, Re_coolant(i));
        Nu_l(i) = ((f / 8) * (Re_coolant(i) - 1000) * Pr_coolant(i)) / (1 + 12.7 * ((f / 8) ^ 0.5) * (Pr_coolant(i) ^ (2/3) - 1)); % Gnielinksy correlation nusselt number [N/A] - 0.5 < Pr < 2000, 3000 < Re < 5e6
        h_l(i) = (Nu_l(i) * k_coolant(i)) / hydraulic_D; % liquid film coefficient [W/m^2-K] (Heister EQ 6.19)

       % land/fin thickness at base
        fin_w(i) = (2*pi*(r_interpolated(i) + t_w_x(i)) - w_c_x(i) * num_channels) / num_channels;
        
        if fin_w(i) <= 0
            error("Invalid channel geometry: fin_w <= 0 at station %d", i);
        end
        
        % Fins: rectangular land/rib, adiabatic tip
        A_c_fin = fin_w(i) * ds_vec(i);      % conduction area
        P_fin = 2 * ds_vec(i);               % two coolant wetted side faces
        
        m_fin = sqrt(h_l(i) * P_fin / (k_w_current(i) * A_c_fin));
        eta_fin(i) = tanh(m_fin * h_c_x(i)) / (m_fin * h_c_x(i));

        A_fin_surf = 2 * h_c_x(i) * ds_vec(i);  % side area only
        A_base = w_c_x(i) * ds_vec(i);       % channel floor area
        
        Qdot_fin(i) = eta_fin(i) * h_l(i) * A_fin_surf * (T_wl(i) - T_coolant(i));
        Qdot_l(i) = h_l(i) * A_base * (T_wl(i) - T_coolant(i)) + Qdot_fin(i);

        % Heat Loop Convergence Check
        if abs(Qdot_g(i) - Qdot_l(i)) > qdot_tolerance && counter < 500

            if (Qdot_g(i) - Qdot_l(i)) > 0
                T_wg_min = T_wg(i);
            else 
                T_wg_max = T_wg(i);
            end 

            T_wg(i) = (T_wg_max + T_wg_min) / 2; 
            counter = counter + 1;
            
            if debug % Debug Tool 
                disp(i)
                disp(T_wg(i))
            end
        else
            Qdot_avg = (Qdot_g(i) + Qdot_l(i)) / 2;
            heatflux(i) = Qdot_avg / A_hot;
            heatflux_fin(i) = Qdot_fin(i) / A_fin_surf;
            
            if debug % Debug Tool 
                Bi_fin(i) = h_l(i) * (A_c_fin / P_fin) / k_w_current(i); % Fin Biot number (should be << 1)
                fin_heat_fraction(i) = Qdot_fin(i) / Qdot_l(i); 
                
                error(i) = Qdot_g(i) - Qdot_l(i);
              
                if coolant_direction
                    if i == 1
                        Qdot_check(i) = mdot_coolant * cp_coolant(i) * (T_coolant(i) - T_amb);
                    else 
                        Qdot_check(i) = mdot_coolant * cp_coolant(i) * (T_coolant(i) - T_coolant(i-1));
                    end
                else
                    if i == steps 
                        Qdot_check(i) = mdot_coolant * cp_coolant(i) * (T_coolant(i) - T_amb);
                    else
                        Qdot_check(i) = mdot_coolant * cp_coolant(i) * (T_coolant(i) - T_coolant(i+1));
                    end
                end
            end

            if coolant_direction
                j = i + 1;
            else
                j = i - 1;
            end
            
            if j >= 1 && j <= steps
                % Update next coolant temperature first
                T_coolant(j) = T_coolant(i) + Qdot_avg / (mdot_channel * cp_coolant(i));
            
                % Carry wall guesses forward
                T_wg(j) = T_wg(i);
                T_wl(j) = T_wl(i);
            
                % Current station geometry
                A_i = h_c_x(i) * w_c_x(i);
                Pwet_i = 2 * h_c_x(i) + 2 * w_c_x(i);
                Dh_i = 4 * A_i / Pwet_i;
            
                % Next station geometry
                A_j = h_c_x(j) * w_c_x(j);
                Pwet_j = 2 * h_c_x(j) + 2 * w_c_x(j);
                Dh_j = 4 * A_j / Pwet_j;
            
                % Estimate next density using updated temperature and current pressure
                if coolant == "isopropyl alcohol"
                    rho_j_est = Data.F_rho(T_K, P_kPa);
                else
                    rho_j_est = double(py.CoolProp.CoolProp.PropsSI('D', 'T', T_coolant(j), 'P', P_coolant(i), coolant));
                end
            
                % Velocities at current and next station
                V_i = mdot_channel / (rho_coolant(i) * A_i);
                V_j = mdot_channel / (rho_j_est * A_j);
            
                % Major loss using channel arc length instead of axial dx
                deltaP_major = PressureLoss_factor * f * rho_coolant(i) / 2 * vel_coolant(i)^2 / hydraulic_D * ds_vec(i);
            
                % Local minor losses
                K_local = 0;
                if i == idx_Cv
                    K_local = K_local + P_minor_coefs(1);
                end
                if i == idx_throat
                    K_local = K_local + P_minor_coefs(2);
                end
                if i == idx_inj
                    K_local = K_local + P_minor_coefs(3);
                end
            
                deltaP_minor = K_local * (rho_coolant(i) * V_i^2 / 2);
            
                % Static pressure change due to velocity and density change
                deltaP_accel = 0.5 * (rho_j_est * V_j^2 - rho_coolant(i) * V_i^2);
            
                % Total pressure update
                P_coolant(j) = P_coolant(i) - deltaP_major - deltaP_minor - deltaP_accel;
            end
                
            Qtot = Qtot + Qdot_avg * num_channels; %total heat through entire chamber [W]
   
            % structural calculations
            if i <= steps 
                yield(i) = interp1(yield_strength(:,1), yield_strength(:,2), T_wg(i), 'linear', 'extrap');
                E_current(i) = interp1(E(:,1), E(:,2), T_wg(i), 'linear', 'extrap');
                CTE_current(i) = interp1(CTE(:,1), CTE(:,2), T_wg(i), 'linear', 'extrap');
                CTE_liq_side(i) = interp1(CTE(:,1), CTE(:,2), T_wl(i), 'linear', 'extrap');
                elong(i) = interp1(elongation_break(:,1), elongation_break(:,2), T_wg(i),'linear','extrap');
                if elong(i) > .25
                    elong(i) = .25;
                end
                epsilon_emax(i) = ((yield(i)*1000000)/ E_current(i));

                deltaT1(i) = T_wg(i) - T_wl(i);
                deltaT2(i) = ((T_wg(i) + T_wl(i))/2) - T_coolant(i); 

                sigma_tp(i) = ( ((P_coolant(i)-P_g(i))/2).*((w_c_x(i)./t_w_x(i)).^2) );
                sigma_tp_cold(i) =  ( ((P_coolant(i))/2).*((w_c_x(i)./t_w_x(i)).^2) );
                sigma_tt(i) = (E_current(i)*CTE_current(i)*heatflux(i)*t_w_x(i))/(2*(1-v)*k_w_current(i));
                sigma_t(i) = ( ((P_coolant(i)-P_g(i))/2).*((w_c_x(i)./t_w_x(i)).^2) ) + (E_current(i)*CTE_current(i)*heatflux(i)*t_w_x(i))/(2*(1-v)*k_w_current(i)); % tangential stress
                sigma_lc(i) = E_current(i)*(CTE_liq_side(i)*(T_wl(i)-T_coolant(i)) + ((CTE_liq_side(i)*deltaT1(i))/(2*(1-v))));
                sigma_ll(i) = E_current(i)*(CTE_current(i)*(T_wg(i)-T_coolant(i)) + ((CTE_current(i)*deltaT1(i))/(2*(1-v))));


                sigma_vc(i) = sqrt(sigma_lc(i)^2 + sigma_t(i)^2 - sigma_lc(i)*sigma_t(i));
                sigma_vl(i) = sqrt(sigma_ll(i)^2 + sigma_t(i)^2 - sigma_ll(i)*sigma_t(i));

                % Calculate total Strains
                epsilon_lc(i) = ((CTE_liq_side(i)*deltaT1(i))/(2*(1-v))) + CTE_liq_side(i)*(T_wl(i)-T_coolant(i));
                epsilon_ll(i) = ((CTE_current(i)*deltaT1(i))/(2*(1-v))) + CTE_current(i)*(T_wg(i)-T_coolant(i));  
                epsilon_tp(i) = ( ((P_coolant(i)-P_g(i))/2).*((w_c_x(i)./t_w_x(i)).^2) )  /E_current(i);
                epsilon_tt(i) = ((E_current(i)*CTE_current(i)*heatflux(i)*t_w_x(i))/(2*(1-v)*k_w_current(i)))  /E_current(i);
                epsilon_t(i) = 1.15 * (( ((P_coolant(i)-P_g(i))/2).*((w_c_x(i)./t_w_x(i)).^2) ) + (E_current(i)*CTE_current(i)*heatflux(i)*t_w_x(i))/(2*(1-v)*k_w_current(i))) / E_current(i); % tangential stress
                epsilon_vc(i) = sqrt(epsilon_lc(i)^2 + epsilon_t(i)^2 - epsilon_lc(i)*epsilon_t(i));
                epsilon_vl(i) = sqrt(epsilon_ll(i)^2 + epsilon_t(i)^2 - epsilon_ll(i)*epsilon_t(i));

                epsilon_tota(i) = ((CTE_current(i)*deltaT1(i))/(2*(1-v))) + CTE_current(i) * deltaT2(i); 
                epsilon_tott(i) = epsilon_t(i);
                epsilon_toteff(i) = (2/sqrt(3)) * sqrt(((epsilon_tott(i)^2)+ epsilon_tott(i)*epsilon_tota(i) + (epsilon_tota(i))^2));
                sigma_a(i) = E_current(i) * epsilon_tota(i);
                sigma_t2(i) = E_current(i) * epsilon_tott(i);
                epsilon_pa(i) = epsilon_tota(i) - epsilon_emax(i);
                epsilon_pt(i) = epsilon_tott(i) - epsilon_emax(i);
                epsilon_emaxeff(i) = (2/sqrt(3)) * sqrt(((epsilon_emax(i)^2) + epsilon_emax(i)*epsilon_emax(i) + (epsilon_emax(i))^2));
                if epsilon_pa(i) < 0
                    epsilon_pa(i) = 0;
                end
                if epsilon_pt(i) < 0
                    epsilon_pt(i) = 0;
                end
                if epsilon_pa(i) == 0 && epsilon_pt(i) == 0 
                    epsilon_eeff(i) = epsilon_toteff(i);
                elseif epsilon_pa(i) > 0 && epsilon_pt(i) ==  0
                    epsilon_eeff(i) = (2/sqrt(3)) * sqrt(((epsilon_tott(i)^2) + epsilon_tott(i)*epsilon_emax(i) + (epsilon_emax(i))^2));
                elseif epsilon_pa(i) == 0 && epsilon_pt(i) > 0 
                    epsilon_eeff(i) = (2/sqrt(3)) * sqrt(((epsilon_emax(i)^2) + epsilon_emax(i)*epsilon_tota(i) + (epsilon_tota(i))^2));
                else 
                    epsilon_eeff(i) = epsilon_emaxeff(i);
                end

                epsilon_peff(i) = (2/sqrt(3)) * sqrt(((epsilon_pt(i)^2) + epsilon_pt(i)*epsilon_pa(i) + (epsilon_pa(i))^2));
                epsilon_cs(i) = 2*yield(i)/E_current(i) + ((elong(i)/2)*((N)^(-1/2)));
                %epsilon_cs_tot(i) = epsilon_cs(i) + epsilon_emax(i); % Total strain after N # of hotfires

                epsilon_cs_tot(i) = epsilon_cs(i) + epsilon_eeff(i)/2; % Total strain after N # of hotfires
                MS_lowcycle(i) = epsilon_cs(i) / (epsilon_peff(i));
                MS(i) = epsilon_cs_tot(i) / (epsilon_toteff(i));
                num_fires(i) = 1/SF * ((2/elong(i))* (epsilon_toteff(i)- (epsilon_eeff(i)/2)))^(-2); %divide by 4 per NASA 5012C
                num_fires_lowcycle(i) = 1/SF * ((2*epsilon_peff(i))/elong(i))^(-2); %divide by 4 per NASA 5012C
                epsilon_cs_spacex(i) = (3.5 * yield(i)*((N)^(-.12)))/E_current(i) + (elong(i)/N)^(.6);
                MS_spacex(i) = epsilon_cs_spacex(i) / epsilon_toteff(i);

                sigma_eff(i) = E_current(i) * epsilon_toteff(i);
                sigma_a(i) = E_current(i) * epsilon_tota(i);
                sigma_t2(i) = E_current(i) * epsilon_tott(i);

            end
                
            converged = 1;
        end

    end % Heat transfer convergence loop end
end % Axial station loop end

overall_MS = min(MS);
overall_MS_lowcycle = min(MS_lowcycle);
overall_MS_spacex = min(MS_spacex);
Engine_life = min(num_fires);
Engine_life_lowcycle = min(num_fires_lowcycle);
yield_SF = min(yield)/(max(sigma_vl)*.000001)-1;
chamber_CDA = mdot_coolant / sqrt(2 * mean(rho_coolant) * (max(P_coolant) - min(P_coolant)));

%% Output
Lifespan = Engine_life;
PressDrop = (max(P_coolant) - min(P_coolant)) / 6894.76;
if DisplayMode == 1
    %% FORMATTED OUTPUT
    fprintf("\nEngine Throttle: %.1f", throttle * 100)
    fprintf("\n\nMargin of safety for engine life of %0.0f hot fires: %.02f", N/8, overall_MS)
    fprintf("\nEngine life (hot fires): %.02f", Engine_life)
    fprintf("\nLowcycle Margin of safety for engine life of %0.0f hot fires: %.02f", N/8, overall_MS_lowcycle)
    fprintf("\nEngine life (hot fires) Lowcycle: %.02f", Engine_life_lowcycle)
    fprintf("\nManson Universal Slopes Margin of safety for engine life of %0.0f hot fires: %.02f", N/8, overall_MS_spacex)
    fprintf("\nSafety factor to yield: %.02f", yield_SF)
    fprintf("\nTotal Heat Input: %.2f kW", Qtot / 10^3)
    fprintf("\nCoolant Pressure Drop: %.2f psi", (max(P_coolant) - min(P_coolant)) / 6894.76)
    fprintf("\nMax Hotwall Temp: %.2f K", max(T_wg))
    fprintf("\nCoolant Exit Temp: %.2f K", max(T_coolant))
    fprintf("\nCoolant Temp Rise: %.2f K", max(T_coolant) - min(T_coolant))
    fprintf("\nChamber CdA: %.2f*10^-5 m^2", chamber_CDA * 10^5)
    
    %% FEA INPUTS
    if FEA_outputs
        ambient_chamber = [x_interpolated; T_g]; % Gas temp [K]
        ambient_coolant = [x_interpolated; T_coolant]; % Coolant temp [K]
        gas_h = [x_interpolated; h_g]; % Hotwall heat transfer coefficient [W/m^2K]
        liquid_h = [x_interpolated; h_l]; % Coldwall heat transfer coefficent [W/m^2K]
        gas_p = [x_interpolated; P_g]; % Gas-side pressure [Pa]
        liquid_p = [x_interpolated; P_coolant]; % Coolant-side pressure [Pa]
        FEA_inputs = [x_interpolated; T_g; T_r; h_g; P_g; T_coolant; h_l; P_coolant];
    
        x_interpolated_fea = -1 * x_interpolated; % fix coordinate system for FEA set up (regen slice FEA)
        Excel_inputs = [h_c_x; w_c_x; x_interpolated_fea; T_g; T_r; h_g; P_g; T_coolant; h_l; P_coolant]';
        
        delete('FEA_regen_large.xls');
        writematrix(Excel_inputs, 'FEA_regen_large.xls');
        % EXCEL COLS: h_c, w_c, x-val, gas-temp, recovery-temp, gas-film-coeff, gas-pressure, coolant-temp, coolant-film-coeff, coolant-pressure 
    end
    
    %% PLOTS
    % Debug Plots
    if debug
        % Heat Flux Error Print
        fprintf("\nHEAT FLUX ERROR %.3f W", sum(error))
        fprintf("\nCOOLANT HEAT INPUT %.3f kW", sum(Qdot_check) / 1000)
        
        % CEA Checker 
        figure()
        plot(x_interpolated, M / sum(M))
        hold on; grid on
        plot(x_interpolated, P_g / sum(P_g))
        plot(x_interpolated, T_g / sum(T_g))
        title("CEA Check Plots")
    
        % Heat Flux Convergence Error
        figure()
        plot(x_interpolated, error)
        grid on
        yline(qdot_tolerance, '-r')
        yline(-qdot_tolerance, '-r')
        title("Heat Loop Convergence Error")
    end
    
    % Temperature Plots
    figure('Name', 'Temperature Plot');
    hold on;
    set(gca, 'FontName', 'Times New Roman')
    yyaxis left
    plot(x_interpolated / 0.0254, T_wg, 'red', 'LineStyle', '-');
    plot(x_interpolated / 0.0254, T_wl, 'magenta', 'LineStyle', '-');
    plot(x_interpolated / 0.0254, T_coolant, 'blue', 'LineStyle', '-')
    plot(x_interpolated / 0.0254, boilingPt, 'cyan', 'LineStyle', '--');
    ylabel('Temperature [K]')
    set(gca, 'YColor', [0 0 0])
    grid on
    yyaxis right
    plot(x_interpolated / 0.0254, r_interpolated / 0.0254, 'green', 'LineStyle', '-');
    ylabel('Radius [in]')
    set(gca, 'YColor', [0 0 0])
    axis auto;
    legend('Gas Wall Temp', 'Liquid Wall Temp','Coolant Temp', 'Boiling Point', 'Contour', 'Location', 'southoutside', 'Orientation', 'horizontal','Location' , 'south')
    title('Wall Temperature Distribution')
    xlabel('Location [in]')
    
    % Coolant Plots
    figure('Name', 'Coolant Properties')
    subplot(2,2,[1,2])
    hold on 
    set(gca, 'FontName', 'Times New Roman')
    yyaxis left
    plot(x_interpolated, P_coolant / 6894.76)
    plot(x_interpolated, P_g / 6894.76, 'y')
    title("Liquid Pressure Loss")
    xlabel("Location [in]")
    ylabel("Pressure [PSI]")
    set(gca, 'YColor', [0 0 0])
    yyaxis right
    plot(x_interpolated, r_interpolated .* u.M2IN, 'green', 'LineStyle', '-');
    ylabel('Radius [in]')
    set(gca, 'YColor', [0 0 0])
    legend("Liquid Pressure","Gas Pressure","Chamber Contour", 'Location', 'southwest')
    grid on
    subplot(2,2,3)
    hold on 
    set(gca, 'FontName', 'Times New Roman')
    plot(x_interpolated, vel_coolant);
    title("Coolant Velocity")
    xlabel("Location [in]")
    ylabel('Velocity [m/s]')
    ylim([0,20])
    grid on
    subplot(2,2,4)
    hold on 
    set(gca, 'FontName', 'Times New Roman')
    plot(x_interpolated, T_coolant);
    title("Coolant Temperature")
    xlabel("Location [in]")
    ylabel('Temperature [K]')
    ylim([275,420])
    grid on
    sgtitle('Coolant Property Plots');
    
    
    % Heat Transfer Plot
    figure('Name', 'Heat Transfer Plots');
    subplot(2,2,[1,2])
    hold on;
    set(gca, 'FontName', 'Times New Roman')
    yyaxis left
    plot(x_interpolated, heatflux / 1000, 'red', 'LineStyle', '-');
    plot(x_interpolated, heatflux_fin / 1000,'yellow', 'LineStyle', '-');
    ylabel('Heat Flux [kW/m^2]')
    set(gca, 'YColor', [0 0 0])
    grid on
    yyaxis right
    plot(x_interpolated, r_interpolated .* u.M2IN, 'green', 'LineStyle', '-');
    ylabel('Radius [in]')
    set(gca, 'YColor', [0 0 0])
    legend('Convective Heat Flux', 'Fin Heat Flux' , 'Chamber Contour','Location','west')
    title('Heat Flux Distribution')
    xlabel('Location [in]')
    subplot(2,2,3)
    hold on 
    set(gca, 'FontName', 'Times New Roman')
    plot(x_interpolated, h_g / 1000)
    title("Gas Heat Transfer Coefficient")
    xlabel("Location [in]");
    ylabel("Bartz HTC [kW/m^2-K]")
    grid on
    subplot(2,2,4)
    hold on 
    set(gca, 'FontName', 'Times New Roman')
    plot(x_interpolated, h_l / 1000)
    title("Liquid Heat Transfer Coefficient")
    xlabel("Location [in]");
    ylabel("Gnielinski HTC [kW/m^2-K]")
    grid on
    
    % Channel Geometry Plots
    figure('Name','Channel Geometry');
    subplot(2,3,[1,2]);
    hold on
    set(gca, 'FontName', 'Times New Roman')
    yyaxis left
    plot(x_interpolated, w_c_x * 1000);
    plot(x_interpolated, h_c_x * 1000);
    title("Channel Dimensions");
    xlabel("Location [in]");
    ylabel("Channel Dimensions [mm]")
    set(gca, 'YColor', [0 0 0])
    yyaxis right
    plot(x_interpolated, r_interpolated .* u.M2IN, 'green', 'LineStyle', '-');
    ylabel('Chamber Contour [in]')
    set(gca, 'YColor', [0 0 0])
    legend('Channel Width', 'Channel Height', 'Chamber Contour','Location','west')
    grid on
    subplot(2,3,3);
    hold on
    set(gca, 'FontName', 'Times New Roman')
    set(gca, 'YAxisLocation', 'right');
    plot(x_interpolated, h_c_x ./ w_c_x);
    title("Channel Aspect Ratio");
    xlabel("Location [in]");
    grid on
    subplot(2,3,[4,5]);
    hold on
    set(gca, 'FontName', 'Times New Roman')
    yyaxis left
    plot(x_interpolated, fin_w * 1000);
    plot(x_interpolated, h_c_x * 1000);
    title("Fin Dimensions")
    xlabel("Location [in]");
    ylabel("Fin Dimensions [mm]")
    set(gca, 'YColor', [0 0 0])
    yyaxis right
    plot(x_interpolated, r_interpolated .* u.M2IN, 'green', 'LineStyle', '-');
    ylabel('Chamber Contour [in]')
    set(gca, 'YColor', [0 0 0])
    legend('Fin Width', 'Fin Height', 'Chamber Contour','Location','west')
    grid on
    subplot(2,3,6);
    hold on
    set(gca, 'FontName', 'Times New Roman')
    set(gca, 'YAxisLocation', 'right');
    plot(x_interpolated, h_c_x ./ fin_w);
    title("Fin Aspect Ratio");
    xlabel("Location [in]");
    grid on
    sgtitle('Channel Geometery Plots');
    
    % Strain Checks
    figure('Name'  , 'Strain Check');
    subplot(2,1,1); hold on; set(gca, 'FontName', 'Times New Roman')
    yyaxis left                               
    plot(x_interpolated, epsilon_toteff * 100);
    plot(x_interpolated, epsilon_cs_tot * 100);
    xlabel("Location [in]");
    ylabel("Strain [%]");
    set(gca, 'YColor', [0 0 0])
    yyaxis right;
    plot(x_interpolated, r_interpolated .* u.M2IN, 'green', 'LineStyle', '-');
    ylabel('Radius [in]');
    set(gca, 'YColor', [0 0 0])
    legend("Total Effective Strain", "Total Allowable Cyclic Strain" , 'Chamber Contour', 'Location' , 'west');
    title('Effective Cyclic Strains')
    grid on
    subplot(2,1,2)
    hold on
    set(gca, 'FontName', 'Times New Roman')
    plot(x_interpolated, epsilon_peff * 100);
    plot(x_interpolated, epsilon_cs * 100) 
    xlabel("Location [in]");
    ylabel("Strain [%]");
    legend("Effective Plastic Strain", "Allowable Cyclic Plastic Strain", 'Location' , 'west')
    title('Plastic Strains')
    grid on
    sgtitle('Cyclic Strains')
    
    % Strain Plots
    figure('Name', 'Structural Strain')
    subplot(2,2,[1,2])
    hold on
    set(gca, 'FontName', 'Times New Roman')
    yyaxis left
    plot(x_interpolated, epsilon_vl * 100);
    plot(x_interpolated, epsilon_vc * 100);
    title("Von Mises Strain")
    xlabel("Location [in]")
    ylabel("Percent Strain [%]")
    set(gca, 'YColor', [0 0 0])
    yyaxis right
    plot(x_interpolated, r_interpolated .* u.M2IN, 'green', 'LineStyle', '-');
    ylabel('Radius [in]')
    set(gca, 'YColor', [0 0 0])
    legend('At the Lands','At the Channels', 'Chamber Contour','Location','west')
    grid on
    subplot(2,2,3)
    hold on 
    set(gca, 'FontName', 'Times New Roman')
    hold on
    plot(x_interpolated, epsilon_t * 100,"b");
    plot(x_interpolated, epsilon_tp * 100,"m");
    plot(x_interpolated, epsilon_tt * 100,"r");
    legend("Total Strain", "Pressure contribution", "Thermal Contribution",'Location','west')
    title("Tangential Strain")
    xlabel("Location [in]")
    ylabel("Percent Tangential Strain [%]")
    hold off
    grid on
    subplot(2,2,4)
    hold on 
    set(gca, 'FontName', 'Times New Roman')
    plot(x_interpolated, epsilon_lc * 100);
    plot(x_interpolated, epsilon_ll * 100);
    title("Longitudinal Strain")
    legend("At the Channels", "At the Lands",'Location','west');
    xlabel("Location [in]")
    ylabel("Percent Longitudinal Strain [%]")
    grid on
    sgtitle('Strain Components')
    
    % Structural Calcs Plot
    figure('Name', 'Structural Stress')
    subplot(2,2,[1,2])
    hold on
    set(gca, 'FontName', 'Times New Roman')
    yyaxis left
    plot(x_interpolated, sigma_vl * 0.000001); 
    plot(x_interpolated, sigma_vc * 0.000001);
    plot(x_interpolated, yield, "yellow");
    title("Von Mises Stress")
    xlabel("Location [in]")
    ylabel("Stress [MPA]")
    set(gca, 'YColor', [0 0 0])
    yyaxis right
    plot(x_interpolated, r_interpolated .* u.M2IN, 'green', 'LineStyle', '-');
    ylabel('Radius [in]')
    set(gca, 'YColor', [0 0 0])
    legend('At the Lands','At the Channels' , 'Yield Stress at Temperature' ,'Chamber Contour','Location','west')
    grid on
    subplot(2,2,3)
    hold on 
    set(gca, 'FontName', 'Times New Roman')
    hold on
    plot(x_interpolated, sigma_t * 0.000001,"b");
    plot(x_interpolated, sigma_tp * 0.000001,"m");
    plot(x_interpolated, sigma_tt * 0.000001,"r");
    plot(x_interpolated, yield, "yellow");
    legend("Total Stress", "Pressure contribution", "Thermal Contribution", 'Yield Stress at Temperature', 'Location','west')
    title("Tangential Stress")
    xlabel("Location [in]")
    ylabel("Stress [MPA]")
    hold off
    grid on
    subplot(2,2,4)
    hold on 
    set(gca, 'FontName', 'Times New Roman')
    plot(x_interpolated, sigma_lc * 0.000001);
    plot(x_interpolated, sigma_ll * 0.000001);
    plot(x_interpolated, yield, "yellow");
    title("Longitudinal Stress")
    legend("At the Channels", "At the Lands", 'Yield Stress at Temperature', 'Location','west');
    xlabel("Location [in]")
    ylabel("Stress [MPA]")
    grid on
    sgtitle('Stress');
    
    % Channel Cross Section Plot
    figure('Name', 'Channel Cross Sections')
    if traditional
        h_c(3) = h_c(2);
        w_c(3) = w_c(2);
        Channel_CS_Plot(num_channels,t_w,h_c,w_c,[r_interpolated(1), r_t, r_interpolated(end)],8,[" Chamber", " Throat", " Nozzle"]);
    else
        Channel_CS_Plot(num_channels,t_w,h_c,w_c,[r_interpolated(1), r_t, r_interpolated(end)],8,[" Chamber", " Throat", " Nozzle"]);
    end
    
    % End De la Code 
    elapsedTime = toc;
    fprintf('\nThe code ran in %.4f seconds.\n', elapsedTime)
end
end

