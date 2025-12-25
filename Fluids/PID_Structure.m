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
[Y0_OX, Y0_FU] = InitializeTanks(TankVolOX, TankVolFU);

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
    Node(2) = Nodes('Tank OX', 2, 0, P_atm, TankVolOX, Y0_OX, 1, 3, false);
    Node(3) = Nodes('Tank IPA', 3, 0, P_atm, TankVolFU, Y0_FU, 2, 4, false);
% Lines
    Node(4) = Nodes('Pre OX', 4, 0, P_atm, 0.00025, [1 0 0], 3, 5, false);
    Node(5) = Nodes('Pre FU', 5, 0, P_atm, 0.00025, [0 1 0], 4, 6, false);
    Node(6) = Nodes('Post OX', 6, 0, P_atm, 0.00025, [0 0 1], 5, 7, false);
    Node(7) = Nodes('Post FU', 7, 0, P_atm, 0.00025, [0 0 1], 6, 8, false);
% Manifold
    Node(8) = Nodes('OX Manifold', 8, 0, P_atm, 0.00025, [0 0 1], 7, 9, false);
    Node(9) = Nodes('FU Manifold', 9, 0, P_atm, 0.00011, [0 0 1], 8, 10, false);
% TADPOLE
    Node(10) = Nodes('TADPOLE', 10, 0, P_atm, 0.00125, [0 0 1], [9 10], 11, true);
% Atmosphere
    Node(11) = Nodes('Atmosphere', 11, 1, P_atm, 1, [0 0 1], 11, [], false);

%% Links
% Tank Press
    Link(1) = ValveLink('Press_OX', 'Throttle', 1, 1, 2, 2.0);
    Link(2) = ValveLink('Press_FU', 'Throttle', 2, 1, 3, 2.0);
% Run Valves
    Link(5) = ValveLink('Main_OX', 'Throttle', 5, 4, 6, 2.9);
    Link(6) = ValveLink('Main_FU', 'Throttle', 6, 5, 7, 2.9);
% Pipe Sections 
    Link(3) = PipeLink('OX Line 1', 3, 2, 4, 0.5, 1e-4, .5);
    Link(4) = PipeLink('FU Line 1', 4, 3, 5, 0.5, 1e-4, .5);
    Link(7) = PipeLink('OX Line 2', 7, 6, 8, 0.5, 1e-4, .5);
    Link(8) = PipeLink('FU Line 2', 8, 7, 9, 0.5, 1e-4, .5);
% Injectors
    Link(9) = ValveLink('Inj OX', 'Orifice', 9, 8, 10, 0.45, 3.146e-5);
    Link(10) = ValveLink('Inj FU', 'Orifice', 10, 9, 10, 0.77, 2.6e-5);
% Nozzle 
    Link(11) = ValveLink('Nozzle', 'Orifice', 11, 10, 11, 0.79, 0.0011);

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



