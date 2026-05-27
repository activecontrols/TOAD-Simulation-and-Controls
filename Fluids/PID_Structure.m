% This file creates the data structure of nodes and links for the lumped
% paramter model of TOAD's P&ID Simulation
% 
%
% Pablo Plata   -   12/02/25


%% System constants and setup
clear;
addpath(".\Fluids\Links and Nodes\")
P_atm = 14.7;
TankVolOX = 0.03075;
TankVolFU = 0.03746;

% Initialize Tank Mass Fractions to 10% Ullage
StartP = 500; %P_atm;
[Y0_OX, Y0_FU] = InitializeTanks(TankVolOX, TankVolFU, StartP);

%% Nodes
% ID Map:
% 1: TK-N2 (COPV)         2: Press Manifold       3: TK-O2-01 (Tank LOX)
% 4: TK-FU-01 (Tank IPA)  5: Pre Main OX          6: Pre Main FU
% 7: Inter OX             8: Inter FU             9: Post Throttle OX
% 10: Post Throttle FU    11: OX Manifold         12: FU Manifold
% 13: SKIPPER (Chamber)   14: Atmosphere          15: Virtual (Igniter Link Node)
% 16: Purge Manifold      17: DART Chamber (New Physical Node)

% Constructor: Nodes(name, ID, Fixed, P0_psi, V_m3, Y0, LinksIN, LinksOUT, isCombustor)
% Source & Pressurization
    Node(1)  = Nodes('TK-N2', 1, 0, 4500, 0.036, [0 0 1], [], 1, false, 'Tank');
    Node(16) = Nodes('Purge Manifold', 16, 0, 550, 1e-3, [0 0 1], 1, [20 15 16], false);
    Node(2)  = Nodes('Press Manifold', 2, 0, P_atm, 1e-3, [0 0 1], 20, [2 3], false);
% Tanks
    Node(3)  = Nodes('TK-O2-01', 3, 0, StartP, TankVolOX, Y0_OX, 2, 5, false, 'Tank');
    Node(4)  = Nodes('TK-FU-01', 4, 0, StartP, TankVolFU, Y0_FU, 2, 6, false, 'Tank');
% Run Lines (Pre-Main)
    Node(5)  = Nodes('Pre Main OX', 5, 0, StartP, 1e-5, [1 0 0], 3, 7, false);
    Node(6)  = Nodes('Pre Main FU', 6, 0, StartP, 1e-5, [0 1 0], 4, 8, false);
% Chill-in
    Node(3).Temp = 90;
    Node(5).Temp = 90;
% Inter Lines (Split flow paths directly into Throttles and DART Solenoids)
    Node(7)  = Nodes('Inter OX', 7, 0, P_atm, 3e-5, [0 0 1], 5, [9 17], false);
    Node(8)  = Nodes('Inter FU', 8, 0, P_atm, 3e-5, [0 0 1], 6, [10 18], false);
% Post-Throttle Lines
    Node(9)  = Nodes('Post Throttle OX', 9, 0, P_atm, 1e-5, [0 0 1], 7, 11, false);
    Node(10) = Nodes('Post Throttle FU', 10, 0, P_atm, 1e-5, [0 0 1], 8, 12, false);
% Manifolds
    Node(11) = Nodes('OX Manifold', 11, 0, P_atm, 1e-5, [0 0 1], 9, 13, false);
    Node(12) = Nodes('FU Manifold', 12, 0, P_atm, 2e-5, [0 0 1], 10, 13, false);
% TADPOLE Main Chamber & DART Torch Chamber
    Node(13) = Nodes('SKIPPER', 13, 0, P_atm, 0.00121, [0 0 1], [11 12 21], 14, true);
    Node(17) = Nodes('DART Chamber', 17, 0, P_atm, 2.7e-5, [0 0 1], [17 18], 21, true);
    System.Combustion.Tau = 0.006;  
% Atmosphere & Spark Mapping Node
    Node(14) = Nodes('Atmosphere', 14, 1, P_atm, 1, [0 0 1], 13, [], false);
    Node(15) = Nodes('Virtual', 15, 1, 0, 1, [0 0 1], 19, [], false);

%% Links
% Pressurization
    Link(1)  = ValveLink('REG-873D', 'Regulator', 1, 1, 16, 0.7, 560, 50, 0.003);
    Link(20) = ValveLink('BV-N2-02', 'Throttle', 20, 16, 2, 2.0);
    Link(2)  = ValveLink('Press OX', 'Check', 2, 2, 3, 25);
    Link(3)  = ValveLink('Press FU', 'Check', 3, 2, 4, 25);
    
% Feed Lines to Mains
    Link(4)  = PipeLink('OX Line 1', 4, 3, 5, 0.4, 5e-4, 20);
    Link(5)  = PipeLink('FU Line 1', 5, 4, 6, 0.4, 5e-4, 20);
    
% Main Valves (03)
    Link(6)  = ValveLink('BV-02-03', 'Solenoid', 6, 5, 7, 2.9);
    Link(7)  = ValveLink('BV-FU-03', 'Solenoid', 7, 6, 8, 2.9);
    
% Throttle Valves (04)
    Link(8)  = ValveLink('BV-02-04', 'Throttle', 8, 7, 9, 2.9);
    Link(9)  = ValveLink('BV-FU-04', 'Throttle', 9, 8, 10, 2.9);
    
% Injector Inertances
    Link(10) = PipeLink('OX Inj Line', 10, 9, 11, 0.01, 5e-4, 20);
    Link(11) = PipeLink('FU Inj Line', 11, 10, 12, 0.01, 5e-4, 20);
    
% Injectors & Nozzle
    Link(12) = ValveLink('Inj OX', 'Orifice', 12, 11, 13, 0.45, 3.146e-5);
    Link(13) = ValveLink('Inj FU', 'Orifice', 13, 12, 13, 0.77, 2.6e-5);
    Link(14) = ValveLink('Nozzle', 'Orifice', 14, 13, 14, 0.79, 0.0011);
    
% N2 Purges
    Link(15) = ValveLink('SV-N2-05', 'Solenoid', 15, 16, 9,  0);  
    Link(16) = ValveLink('SV-N2-06', 'Solenoid', 16, 16, 10, 0);
    
% DART Runs
    Link(17) = ValveLink('SV-DART-OX', 'Solenoid', 17, 7, 17, 0); % Matched to 0.063" Orifice Area
    Link(18) = ValveLink('SV-DART-FU', 'Solenoid', 18, 8, 17, 0); % Matched to 0.029" Orifice Area
    Link(22) = ValveLink('SV-N2-07', 'Solenoid', 22, 16, 17, 0);  % Igniter Purge
    
% Spark Plug Command
    Link(19) = ValveLink('Spark', 'Signal', 19, 15, 15, 0);

% DART Flame Tube Nozzle Link
    Link(21) = ValveLink('DART Nozzle', 'Orifice', 21, 17, 13, 1.0, 2.62e-5);

%% Pre-processing (subdividing into Dynamic and Algebraic Links)
isDynamic = strcmp({Link.Type}, 'Pipe');
DynamicLink = Link(isDynamic);
AlgebraicLink = Link(~isDynamic);

%% Define System
System.Nodes = Node;
System.Links.Dynamic = DynamicLink;
System.Links.Algebraic = AlgebraicLink;
System.Constants.C = 1250;
MALength = size(System.Links.Algebraic, 2);
MDLength = size(System.Links.Dynamic, 2);

% Densities, (OX, FU, GN2) in kg/m^3
System.Constants.RhoArray = [1141, 786, 45];

% Build Link Map
LinkMap = zeros(MALength + MDLength, 2);
for i = 1:MDLength
    L = System.Links.Dynamic(i);
    LinkMap(L.ID, :) = [L.Up, L.Down];
end
for i = 1:MALength
    L = System.Links.Algebraic(i);
    LinkMap(L.ID, :) = [L.Up, L.Down];
end
System.LinkMap = LinkMap;

% Pre-cache Igniter ID and combustor proprieties
System.Combustion.IgniterID = find(strcmp({System.Links.Algebraic.Name}, 'Spark'), 1);
System.Combustion.CEA_OF = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, ...
              1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, ...
              2.8, 2.9, 3.0, 3.1, 3.2, 3.3];

System.Combustion.CEA_Temp = [1562, 1879, 2172, 2436, 2668, 2862, 3017, 3135, 3219, 3277, ...
                3315, 3338, 3352, 3359, 3361, 3360, 3356, 3349, 3341, 3332, ...
                3322, 3311, 3300, 3288, 3275, 3262];

% Enthalpy handling
System.Constants.T_Table = 70:10:400;
T = System.Constants.T_Table;

% N2
t_n2 = T / 1000;
MW_N2 = 0.0280134; % kg/mol
A=28.98641; B=1.853978; C=-9.647459; D=16.63537; E=0.000117; F=-8.671914; H=0.0;
Cp_mol_N2 = A + B*t_n2 + C*t_n2.^2 + D*t_n2.^3 + E./t_n2.^2;
h_mol_N2 = A*t_n2 + B.*t_n2.^2/2 + C.*t_n2.^3/3 + D.*t_n2.^4/4 - E./t_n2 + F - H;

System.Constants.Cp_N2_Table = Cp_mol_N2 / MW_N2;          % J/(kg-K)
System.Constants.h_N2_Table  = (h_mol_N2 * 1000) / MW_N2;  % J/kg

% LOX wrt 70K
MW_OX = 0.031998; % kg/mol
Cp_mol_OX = 53.0 + (T - 70.0) * (6.0 / 50.0);
h_mol_OX  = 53.0 * (T - 70.0) + (0.12 / 2.0) * ((T - 70.0).^2);

System.Constants.Cp_OX_Table = Cp_mol_OX / MW_OX;
System.Constants.h_OX_Table  = h_mol_OX / MW_OX;

% IPA wrt 180K
MW_FU = 0.060096; % kg/mol
Cp_mol_FU = 161.2 + (T - 298.15) * 0.2;
h_mol_FU  = 161.2 * (T - 180.0) + (0.2 / 2.0) * ((T - 298.15).^2 - (180.0 - 298.15).^2);

System.Constants.Cp_FU_Table = Cp_mol_FU / MW_FU;
System.Constants.h_FU_Table  = h_mol_FU / MW_FU;

%% Valve Manager Configuration
Scheduler = ValveManager(System);
Scheduler.SetTau('BV-N2-02', 0.05);
Scheduler.SetTau('BV-02-03', 0.05); 
Scheduler.SetTau('BV-FU-03', 0.05);
Scheduler.SetTau('BV-02-04', 0.07); 
Scheduler.SetTau('BV-FU-04', 0.07);
Scheduler.SetTau('SV-N2-05', 0.05); 
Scheduler.SetTau('SV-N2-06', 0.05); 
Scheduler.SetTau('SV-N2-07', 0.1); 

% Initialize Time Constraints for Torch Run Valves
Scheduler.SetTau('SV-DART-OX', 0.04); 
Scheduler.SetTau('SV-DART-FU', 0.04);
Scheduler.SetTau('Spark', 0.001);

% Test Autosequence Matrix Mapping
[Valves1(1), Valves1(2)] = ValveController(250, 1.2, 550);

AutoSequence = {
    % Pressurize Main Run Lines
    0.0,    'BV-N2-02',       2.0; 

    % Open Mains
    3.0,    'BV-02-03',       2.9;
    3.0,    'BV-FU-03',       2.9;
    
   % Spark
    3.8,    'Spark',       1.0;
    4.1,   'SV-DART-OX',     0.014;
    4.1,   'SV-DART-FU',     0.0195;

    % DART shutdown
    4.50,   'Spark',          0.0;
    5.00,   'SV-DART-OX',     0.0;
    5.00,   'SV-DART-FU',     0.0;
    
    % Main Chamber Startup
    4.2,    'BV-02-04',       Valves1(1);
    4.2,    'BV-FU-04',       Valves1(2);
    
    % Nominal Engine Cutoff Sequence
    6.0,    'BV-02-03',       0.0;
    6.0,    'BV-FU-03',       0.0;
    6.0,    'BV-02-04',       0.0;
    6.0,    'BV-FU-04',       0.0;

    % System Safing Nitrogen Purges
    6.07,   'SV-N2-05',       0.2;
    6.07,   'SV-N2-06',       0.2;
    6.07,   'SV-N2-07',       0.1;
    10.0,   'SV-N2-05',       0.0;
    10.0,   'SV-N2-06',       0.0;
    10.0,   'BV-N2-02',       0.0;
    10.0,   'SV-N2-07',       0.0;
};

Scheduler.LoadSequence(AutoSequence);

