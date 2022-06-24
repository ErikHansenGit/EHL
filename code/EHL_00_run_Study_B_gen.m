close all; clc; clearvars;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EHL study A setup and run code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main developer:
% Erik Hansen: erik.hansen@kit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information:
% This script automatically executes all setups and subsequent simulations
% of study B to replicate the setup of Bertocchi et al., 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [-] amount of different simulations to be run:
autorun.sim.N = 8;      
% Set parameters which are different for the simulations within this study:
% Study and Simualtion Identifiers:
autorun.alg.Sim_ID{1} = 'Study_B/Nr_1/';
autorun.alg.Sim_ID{2} = 'Study_B/Nr_2/';
autorun.alg.Sim_ID{3} = 'Study_B/Nr_3/';
autorun.alg.Sim_ID{4} = 'Study_B/Nr_4/';
autorun.alg.Sim_ID{5} = 'Study_B/Nr_5/';
autorun.alg.Sim_ID{6} = 'Study_B/Nr_6/';
autorun.alg.Sim_ID{7} = 'Study_B/Nr_7/';
autorun.alg.Sim_ID{8} = 'Study_B/Nr_8/';

% [-]   flag whether solid body is elastic on macroscopic scale (true) or rigid (false)
autorun.sld.flag_elastic(1)   = false;     
autorun.sld.flag_elastic(2)   = true; 
autorun.sld.flag_elastic(3)   = false; 
autorun.sld.flag_elastic(4)   = true; 
autorun.sld.flag_elastic(5)   = false;     
autorun.sld.flag_elastic(6)   = true; 
autorun.sld.flag_elastic(7)   = false; 
autorun.sld.flag_elastic(8)   = true; 

% [m]       domain width
autorun.geo.L_x2(1:2)            = 10e-3;   
autorun.geo.L_x2(3:4)            = 30*10e-3; 
autorun.geo.L_x2(5:6)            = 10e-3;   
autorun.geo.L_x2(7:8)            = 30*10e-3; 
% [m]      pocket width
autorun.geo.l_x2(1:2)            = 7.0e-3;                          
autorun.geo.l_x2(3:4)            = 30*7.0e-3;      
autorun.geo.l_x2(5:6)            = 7.0e-3;                          
autorun.geo.l_x2(7:8)            = 30*7.0e-3;

% [-]   flag, which discretization scheme is used for the Couette term 
flag_Couette_Discrete_1 = 1; % 1.Order Upwind
flag_Couette_Discrete_2 = 2; % 2.Order QUICK
autorun.alg.FBNS.flag_Couette_Discrete(1:4) = flag_Couette_Discrete_1;  
autorun.alg.FBNS.flag_Couette_Discrete(5:8) = flag_Couette_Discrete_2;     
    
autorun.sim.it = 1;
while autorun.sim.it <= autorun.sim.N
    run EHL_01_setup_Study_B_gen.m
    run EHL_02_mainprocess.m
    autorun.sim.it = autorun.sim.it + 1;
end