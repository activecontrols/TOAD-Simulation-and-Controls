% This file creates the data structure of nodes and links for the lumped
% paramter model of TOAD's P&ID Simulation
%
% Pablo Plata   -   12/02/25


%% System constants and setup
clear;
addpath(".\Links and Nodes\")
rho_N2 = 1000;
rho_OX = 1000;
rho_FU = 1000;
P_atm = 14.7;
C = 1500;

%% Nodes
% Nodes are defined as the control volumes between diffrent links. Each
% node has a fixed volume V_k. For the purposes of this Simulation, we will
% assign it a length and an area, and constrain the nodes to straight,
% constant area pipes. Each node will also have a defined speed of sound 
% across it, aswell as a friction loss coefficient. Each node
% also has a parameter listing the connecting links. A pressure state is
% assigned to each node.

Node(1) = Nodes('Tank', 1, 0, 550, 0.35, [0 1 0], [], 1);
Node(2) = Nodes('Up', 2, 0, P_atm, 0.001, [0 0 1], 1, 2);
Node(3) = Nodes('Down', 3, 0, P_atm, 0.005, [0 0 1], 2, 3);
Node(4) = Nodes('Atm', 4, 1, P_atm, 0, [0 0 1], 3, []);

%% Links
% Links are defined as the flow elements that connect nodes. There are
% diffrent types of links (ex: "Pipe", "Throttle_Valve", "Binary_Valve",
% "Orifice", "Check_Valve"...). Each link as assigned node indexes it
% connects, aswell as diffrent parameters. If the link is a pipe section,
% it has defined geometric parameters and an assigned ODE. Else, the link
% is assigned an algebraic massflow equation depending on it's type (ex: Cv
% for throttle valves, Orifice equation for Orifices...).

Link(1) = PipeLink('Pipe1', 1, 1, 2, rho_FU, 2, 1e-4, 3);
Link(2) = ValveLink('Valve', 'Throttle', 2, 2, 3, rho_FU, 2);
Link(3) = PipeLink('Pipe2', 3, 3, 4, rho_FU, 0.2, 1e-4, 3);

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
System.Constants.RhoArray = [1141, 786, 1.25];

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



