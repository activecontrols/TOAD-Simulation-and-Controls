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
% 1: Regulator
% 2: Tank LOX
% 3: Tank IPA
% 4: Line LOX (Post-Valve)
% 5: Line IPA (Post-Valve)
% 6: Combustor
% 7: Atmosphere
% Constructor: Nodes(name, ID, Fixed, P0_psi, V_m3, Y0, LinksIN, LinksOUT, isCombustor)
% Tanks
    Node(1) = Nodes('Regulator', 1, 1, 550, 5, [0 0 1], [], [1 2], false);
    Node(2) = Nodes('Tank OX', 2, 0, StartP, TankVolOX, Y0_OX, 1, 3, false, 'Tank');
    Node(3) = Nodes('Tank IPA', 3, 0, StartP, TankVolFU, Y0_FU, 2, 4, false, 'Tank');
% Lines
    Node(4) = Nodes('Pre OX', 4, 0, StartP, 0.0001, [1 0 0], 3, 5, false);
    Node(5) = Nodes('Pre FU', 5, 0, StartP, 0.0001, [0 1 0], 4, 6, false);
    Node(6) = Nodes('Post OX', 6, 0, P_atm, 0.0001, [0 0 1], 5, 7, false);
    Node(7) = Nodes('Post FU', 7, 0, P_atm, 0.0001, [0 0 1], 6, 8, false);
% Manifold
    Node(8) = Nodes('OX Manifold', 8, 0, P_atm, 1e-4, [0 0 1], 7, 9, false);
    Node(9) = Nodes('FU Manifold', 9, 0, P_atm, 1e-4, [0 0 1], 8, 10, false);
% TADPOLE
    Node(10) = Nodes('TADPOLE', 10, 0, P_atm, 0.00121, [0 0 1], [9 10], 11, true);
    System.Combustion.Tau = 0.006;  % Combustion deadtime parameter
% Atmosphere
    Node(11) = Nodes('Atmosphere', 11, 1, P_atm, 1, [0 0 1], 11, [], false);
% Virtual Node (added for igniter connection)
    Node(12) = Nodes('Virtual', 12, 1, 0, 1, [0 0 1], [], [], false);

%% Links
% Tank Press
    Link(1) = ValveLink('Press OX', 'Throttle', 1, 1, 2, 2.0);
    Link(2) = ValveLink('Press FU', 'Throttle', 2, 1, 3, 2.0);
% Run Valves
    Link(5) = ValveLink('Main OX', 'Throttle', 5, 4, 6, 2.9);
    Link(6) = ValveLink('Main FU', 'Throttle', 6, 5, 7, 2.9);
% Pipe Sections 
    Link(3) = PipeLink('OX Line 1', 3, 2, 4, 0.5, 5e-4, 5);
    Link(4) = PipeLink('FU Line 1', 4, 3, 5, 0.5, 5e-4, 5);
% Injector Inertances
    Link(7) = PipeLink('OX Inj 1', 7, 6, 8, 0.3, 5e-4, 5);
    Link(8) = PipeLink('FU Inj 1', 8, 7, 9, 0.3, 5e-4, 5);
% Injectors
    Link(9) = ValveLink('Inj OX', 'Orifice', 9, 8, 10, 0.45, 3.146e-5);
    Link(10) = ValveLink('Inj FU', 'Orifice', 10, 9, 10, 0.77, 2.6e-5);
% Nozzle 
    Link(11) = ValveLink('Nozzle', 'Orifice', 11, 10, 11, 0.79, 0.0011);
% N2 Purges
    Link(12) = ValveLink('Purge OX', 'Check', 12, 1, 6, 0);
    Link(13) = ValveLink('Purge FU', 'Check', 13, 1, 7, 0);
% Igniter Link
    Link(14) = ValveLink('Igniter', 'Signal', 14, 12, 12, 0);

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
System.Combustion.IgniterID = find(strcmp({System.Links.Algebraic.Name}, 'Igniter'), 1);
System.Combustion.CEA_OF = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, ...
              1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, ...
              2.8, 2.9, 3.0, 3.1, 3.2, 3.3];

System.Combustion.CEA_Temp = [1562, 1879, 2172, 2436, 2668, 2862, 3017, 3135, 3219, 3277, ...
                3315, 3338, 3352, 3359, 3361, 3360, 3356, 3349, 3341, 3332, ...
                3322, 3311, 3300, 3288, 3275, 3262];

%% Valve Manager Tests
Scheduler = ValveManager(System);
Scheduler.SetTau('Press OX', 0.005);
Scheduler.SetTau('Press FU', 0.005);
Scheduler.SetTau('Purge OX', 0.005);
Scheduler.SetTau('Purge FU', 0.005);
Scheduler.SetTau('Main OX', 0.07);
Scheduler.SetTau('Main FU', 0.07);
Scheduler.SetTau('Igniter', 0.001);

% Test Autosequence for System
[Valves1(1), Valves1(2)] = ValveController(250, 1.2, 550);
[Valves2(1), Valves2(2)] = ValveController(85, 1.2, 550);
AutoSequence = {
    % Open Press Lines
    0.0,    'Press OX',     0.7;
    0.0,    'Press FU',     0.7;
    
    % Starup 
    6.0,  'Main FU',      Valves1(2);
    6.0,  'Main OX',      Valves1(1);

    % Igniter
    6.03,  'Igniter',      1;

    % Throttle Down
    8.0,    'Main FU',      Valves2(2);
    8.0,    'Main OX',      Valves2(1);

    % Throttle Up
    13.0,   'Main FU',      Valves1(2);
    13.0,   'Main OX',      Valves1(1);
    
    % Shutdown
    16.0,   'Main FU',      0.0;
    16.0,   'Main OX',      0.0;
    16.3,   'Igniter',      0.0;

    % Purges
    16.07,   'Purge OX',     .5;
    16.07,   'Purge FU',     .5;
    
    % Purges and Press shutdown
    17.0,   'Purge OX',     0.0;
    17.0,   'Purge FU',     0.0;
    17.0,   'Press OX',     0.0;
    17.0,   'Press FU',     0.0;
};

%% Purge Test Sequence
% AutoSequence = {
%     % Open Press Lines
%     0.0,    'Press OX',     2.0;
%     0.0,    'Press FU',     2.0;
% 
%     % Purge Test
%     7,   'Purge OX',     .2;
%     7,   'Purge FU',     .2;
%     9,   'Purge OX',      0;
%     9,   'Purge FU',      0;
% };

Scheduler.LoadSequence(AutoSequence);



