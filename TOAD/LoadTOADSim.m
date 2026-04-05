%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file loads in all the constants and parameters for the Simulink into
% workspace. Please always run this file before running a full-scale
% simulation if you've made any changes to trajectory, controls, filtering,
% or others.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize parameters and clear functions
% Initial conditions for state
clear;
clear ref_generator3;
clear inputfcn3;
clear EstimateStateFCN;
clear SensorSimulation;
clear GPS_Sim;
clear DigitalNF;
constants_port;

%% Create constants struct for TOAD (Approximate values, all metric)
constantsTOAD.m_dry = 132.5;
constantsTOAD.g = 9.80145; 
constantsTOAD.rTB = 1.1;
constantsTOAD.J = diag([80 80 15]);
constantsTOAD.MaxThrust = 2446.52;
constantsTOAD.MaxMdot = 1.3204;
constantsTOAD.OF = 1;
constantsTOAD.OxMass = 20.78;   constantsTOAD.FuMass = 19.79;
constantsTOAD.OxHeight = 0.377; constantsTOAD.FuHeight = 0.495;
constantsTOAD.OxRadius = 0.146; constantsTOAD.FuRadius = 0.146;
constantsTOAD.Ox_Z = 0.85;      constantsTOAD.Fu_Z = 1.35;




