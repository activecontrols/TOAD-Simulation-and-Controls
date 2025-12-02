% This file creates the data structure of nodes and links for the lumped
% paramter model of TOAD's P&ID Simulation
%
% Pablo Plata   -   12/02/25

%% Nodes
% Nodes are defined as the control volumes between diffrent links. Each
% node has a fixed volume V_k. For the purposes of this Simulation, we will
% assign it a length and an area, and constrain the nodes to straight,
% constant area pipes. Each node will also have a defined speed of sound 
% across it, aswell as a friction loss coefficient. Each node
% also has a parameter listing the connecting links. A pressure state is
% assigned to each node.

%% Links
% Links are defined as the flow elements that connect nodes. There are
% diffrent types of links (ex: "Pipe", "Throttle_Valve", "Binary_Valve",
% "Orifice", "Check_Valve"...). Each link as assigned node indexes it
% connects, aswell as diffrent parameters. If the link is a pipe section,
% it has defined geometric parameters and an assigned ODE. Else, the link
% is assigned an algebraic massflow equation depending on it's type (ex: Cv
% for throttle valves, Orifice equation for Orifices...).

