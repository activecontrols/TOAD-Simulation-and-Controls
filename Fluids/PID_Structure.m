% This file creates the data structure of nodes and links for the lumped
% paramter model of TOAD's P&ID Simulation
%
% Pablo Plata   -   12/02/25


%% System constants and setup
clear;
addpath(".\Fluids\Links and Nodes\")
rho_FU = 785;
rho_OX = 1140;
P_atm = 14.7;

R_N2_init = 296.8;  % J/kg*K
T_amb_init = 293;   % K
P_init_Pa = P_atm * 6895; % Convert 14.7 psi to Pa for density calc

TargetUllage = 0.1; 
TankVol = 0.0283;   

% Calculate Phase Volumes
Vol_Gas = TankVol * TargetUllage;
Vol_Liq = TankVol * (1 - TargetUllage);

% Calculate Initial Gas Density (Ideal Gas Law)
Rho_Gas_Init = P_init_Pa / (R_N2_init * T_amb_init);

% Calculate Component Masses
Mass_Gas = Vol_Gas * Rho_Gas_Init;
Mass_Liq = Vol_Liq * rho_FU;
Mass_Total = Mass_Gas + Mass_Liq;

% Calculate Mass Fractions Y = [OX, FU, N2]
Y_N2 = Mass_Gas / Mass_Total;
Y_FU = Mass_Liq / Mass_Total;

% Print Check
fprintf('Initializing Tank with %.1f%% Ullage.\n', TargetUllage*100);
fprintf('Mass Fractions -> FU: %.6f, N2: %.6f\n', Y_FU, Y_N2);

%% Nodes
% Nodes are defined as the control volumes between diffrent links. Each
% node has a fixed volume V_k. For the purposes of this Simulation, we will
% assign it a length and an area, and constrain the nodes to straight,
% constant area pipes. Each node will also have a defined speed of sound 
% across it, aswell as a friction loss coefficient. Each node
% also has a parameter listing the connecting links. A pressure state is
% assigned to each node.

Node(1) = Nodes('Regulator', 1, 1, 550, 5, [0 0 1], [], 1, false);
Node(2) = Nodes('Tank', 2, 0, P_atm, TankVol, [0 Y_FU Y_N2], 1, 2, false);
Node(3) = Nodes('Up', 3, 0, P_atm, 0.001, [0 Y_FU Y_N2], 2, 3, false);
Node(4) = Nodes('Down', 4, 0, P_atm, 0.001, [0 0 1], 3, 4, false);
Node(5) = Nodes('Outlet', 5, 1, P_atm, 0, [0 0 1], 4, [], false);

%% Links
% Links are defined as the flow elements that connect nodes. There are
% diffrent types of links (ex: "Pipe", "Throttle_Valve", "Binary_Valve",
% "Orifice", "Check_Valve"...). Each link as assigned node indexes it
% connects, aswell as diffrent parameters. If the link is a pipe section,
% it has defined geometric parameters and an assigned ODE. Else, the link
% is assigned an algebraic massflow equation depending on it's type (ex: Cv
% for throttle valves, Orifice equation for Orifices...).

Link(1) = ValveLink('COPV Valve', 'Throttle', 1, 1, 2, rho_FU, 2);
Link(2) = PipeLink('Tank2Valve', 2, 2, 3, rho_FU, 1, 1e-4, 15);
Link(3) = ValveLink('Outlet Valve', 'Throttle', 3, 3, 4, rho_FU, 0);
Link(4) = PipeLink('Valve2Outlet', 4, 4, 5, rho_FU, 0.4, 1e-4, 15);

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



