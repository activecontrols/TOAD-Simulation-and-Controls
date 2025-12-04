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

Node(1) = Nodes('Tank', 1, 0, 550, [], 1);
Node(2) = Nodes('Up', 2, 0, P_atm, 1, 2);
Node(3) = Nodes('Down', 3, 0, P_atm, 2, 3);
Node(4) = Nodes('Atm', 4, 1, P_atm, 3, []);

%% Links
% Links are defined as the flow elements that connect nodes. There are
% diffrent types of links (ex: "Pipe", "Throttle_Valve", "Binary_Valve",
% "Orifice", "Check_Valve"...). Each link as assigned node indexes it
% connects, aswell as diffrent parameters. If the link is a pipe section,
% it has defined geometric parameters and an assigned ODE. Else, the link
% is assigned an algebraic massflow equation depending on it's type (ex: Cv
% for throttle valves, Orifice equation for Orifices...).

Link(1) = PipeLink('Pipe1', 1, 1, 2, 0.3, 1e-4, 2, rho_FU);
Link(2) = ValveLink('Valve', 'Throttle', 2, 2, 3, rho_FU, 2);
Link(3) = PipeLink('Pipe2', 3, 3, 4, 0.3, 1e-4, 2, rho_FU);

%% Define System
System.Nodes = Node;
System.Links = Link;


