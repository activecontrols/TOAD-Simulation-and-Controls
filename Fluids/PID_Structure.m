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
% 1: TK-N2
% 2: Press Manifold
% 3: TK-02-01 (Tank LOX)
% 4: TK-FU-01 (Tank IPA)
% 5: Pre Main OX 
% 6: Pre Main FU 
% 7: Inter OX (Between Main 03 and Throttle 04)
% 8: Inter FU (Between Main 03 and Throttle 04)
% 9: Post Throt OX 
% 10: Post Throt FU 
% 11: OX Manifold
% 12: FU Manifold
% 13: SKIPPER (Chamber)
% 14: Atmosphere
% 15: DART (Igniter)

% Constructor: Nodes(name, ID, Fixed, P0_psi, V_m3, Y0, LinksIN, LinksOUT, isCombustor)
% Source & Pressurization
    Node(1)  = Nodes('TK-N2', 1, 1, 550, 5, [0 0 1], [], [], false);
    Node(2)  = Nodes('Press Manifold', 2, 0, P_atm, 1e-4, [0 0 1], 1, [3 4], false);
% Tanks
    Node(3)  = Nodes('TK-O2-01', 3, 0, StartP, TankVolOX, Y0_OX, 2, 5, false, 'Tank');
    Node(4)  = Nodes('TK-FU-01', 4, 0, StartP, TankVolFU, Y0_FU, 2, 6, false, 'Tank');
% Run Lines (Pre-Main)
    Node(5)  = Nodes('Pre Main OX', 5, 0, StartP, 1e-5, [1 0 0], 3, 7, false);
    Node(6)  = Nodes('Pre Main FU', 6, 0, StartP, 1e-5, [0 1 0], 4, 8, false);
% Chill-in
    Node(3).Temp = 90;
    Node(5).Temp = 90;
% Inter Lines (Between 03 and 04)
    Node(7)  = Nodes('Inter OX', 7, 0, P_atm, 1e-5, [0 0 1], 5, 9, false);
    Node(8)  = Nodes('Inter FU', 8, 0, P_atm, 1e-5, [0 0 1], 6, 10, false);
% Post-Throttle Lines
    Node(9)  = Nodes('Post Throttle OX', 9, 0, P_atm, 1e-5, [0 0 1], 7, 11, false);
    Node(10) = Nodes('Post Throttle FU', 10, 0, P_atm, 1e-5, [0 0 1], 8, 12, false);
% Manifolds
    Node(11) = Nodes('OX Manifold', 11, 0, P_atm, 2e-5, [0 0 1], 9, 13, false);
    Node(12) = Nodes('FU Manifold', 12, 0, P_atm, 2e-5, [0 0 1], 10, 13, false);
% TADPOLE
    Node(13) = Nodes('SKIPPER', 13, 0, P_atm, 0.00121, [0 0 1], [11 12], 14, true);
    System.Combustion.Tau = 0.006;  
% Atmosphere & Virtual
    Node(14) = Nodes('Atmosphere', 14, 1, P_atm, 1, [0 0 1], 13, [], false);
    Node(15) = Nodes('Virtual', 15, 1, 0, 1, [0 0 1], [], [], false);

%% Links
% Pressurization
    Link(1) = ValveLink('BV-N2-02', 'Throttle', 1, 1, 2, 2.0); 
    Link(2) = PipeLink('Press Line OX', 2, 2, 3, 0.5, 5e-4, 50); % Manifold to Tank
    Link(3) = PipeLink('Press Line FU', 3, 2, 4, 0.5, 5e-4, 50); % Manifold to Tank
    
% Feed Lines to Mains
    Link(4) = PipeLink('OX Line 1', 4, 3, 5, 0.4, 5e-4, 300);
    Link(5) = PipeLink('FU Line 1', 5, 4, 6, 0.4, 5e-4, 300);
    
% Main Valves (03) - Modeled as dynamic Checks
    Link(6) = ValveLink('BV-02-03', 'Check', 6, 5, 7, 2.9);
    Link(7) = ValveLink('BV-FU-03', 'Check', 7, 6, 8, 2.9);
    
% Throttle Valves (04) - Flow control stage
    Link(8) = ValveLink('BV-02-04', 'Throttle', 8, 7, 9, 2.9);
    Link(9) = ValveLink('BV-FU-04', 'Throttle', 9, 8, 10, 2.9);
    
% Injector Inertances
    Link(10) = PipeLink('OX Inj Line', 10, 9, 11, 0.3, 5e-4, 300);
    Link(11) = PipeLink('FU Inj Line', 11, 10, 12, 0.3, 5e-4, 300);
    
% Injectors & Nozzle
    Link(12) = ValveLink('Inj OX', 'Orifice', 12, 11, 13, 0.45, 3.146e-5);
    Link(13) = ValveLink('Inj FU', 'Orifice', 13, 12, 13, 0.77, 2.6e-5);
    Link(14) = ValveLink('Nozzle', 'Orifice', 14, 13, 14, 0.79, 0.0011);
    
% N2 Purges (Straight from Source)
    Link(15) = ValveLink('SV-N2-05', 'Check', 15, 1, 9, 0);  % To Post Throt
    Link(16) = ValveLink('SV-N2-06', 'Check', 16, 1, 10, 0); % To Post Throt
    Link(17) = ValveLink('SV-N2-07', 'Check', 17, 1, 11, 0);  % To Manifold
    Link(18) = ValveLink('SV-N2-08', 'Check', 18, 1, 12, 0);  % To Manifold
    
% Igniter Link
    Link(19) = ValveLink('DART Ign', 'Signal', 19, 15, 15, 0);

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
System.Combustion.IgniterID = find(strcmp({System.Links.Algebraic.Name}, 'DART Ign'), 1);
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

%% Valve Manager Tests
Scheduler = ValveManager(System);
Scheduler.SetTau('BV-N2-02', 0.05);
Scheduler.SetTau('BV-02-03', 0.05); % Mains open fast
Scheduler.SetTau('BV-FU-03', 0.05);
Scheduler.SetTau('BV-02-04', 0.07); % Throttles articulate slightly slower
Scheduler.SetTau('BV-FU-04', 0.07);
Scheduler.SetTau('SV-N2-05', 0.05); % Post OX Purge
Scheduler.SetTau('SV-N2-06', 0.05); % Post FU Purge
Scheduler.SetTau('SV-N2-07', 0.05); % Ign/Man OX Purge
Scheduler.SetTau('SV-N2-08', 0.05); % Ign/Man FU Purge
Scheduler.SetTau('DART Ign', 0.05);

% Test Autosequence for System
[Valves1(1), Valves1(2)] = ValveController(250, 1.2, 550);
[Valves2(1), Valves2(2)] = ValveController(85, 1.2, 550);

AutoSequence = {
    % Pressurize System
    0.0,    'BV-N2-02',       2.0; 

    % Open Throttles to set angle
    3.0,    'BV-02-04',       Valves1(1);
    3.0,    'BV-FU-04',       Valves1(2);
    
    % Startup (Mains Open 100%, Throttles to Target)
    6.0,    'BV-02-03',       2.9;
    6.0,    'BV-FU-03',       2.9;

    % Igniter
    6.03,   'DART Ign',       1.0;

    % Throttle Down (Mains stay open, Throttles adjust)
    8.0,    'BV-02-04',       Valves2(1);
    8.0,    'BV-FU-04',       Valves2(2);

    % Throttle Up
    13.0,   'BV-02-04',       Valves1(1);
    13.0,   'BV-FU-04',       Valves1(2);
    
    % Shutdown (Close both Mains and Throttles)
    16.0,   'BV-02-03',       0.0;
    16.0,   'BV-FU-03',       0.0;
    16.0,   'BV-02-04',       0.0;
    16.0,   'BV-FU-04',       0.0;
    16.3,   'DART Ign',       0.0;

    % All 4 Purges On
    16.07,  'SV-N2-05',       0.5;
    16.07,  'SV-N2-06',       0.5;
    16.07,  'SV-N2-07',       0.5;
    16.07,  'SV-N2-08',       0.5;
    
    % Purges and Press shutdown
    17.0,   'SV-N2-05',       0.0;
    17.0,   'SV-N2-06',       0.0;
    17.0,   'SV-N2-07',       0.0;
    17.0,   'SV-N2-08',       0.0;
    17.0,   'BV-N2-02',       0.0;
};

Scheduler.LoadSequence(AutoSequence);