close all; clc; clearvars;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EHL study A setup and run code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main developer:
% Erik Hansen: erik.hansen@kit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information:
% This script automatically executes all setups and subsequent simulations
% of study A to replicate the pseudo 1d setup of the inclined slider with 
% rectangular pocket.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [-] amount of different simulations to be run:
autorun.sim.N = 20;      
% Set parameters which are different for the simulations within this study:
% Study and Simualtion Identifiers:
autorun.alg.Sim_ID{1}  = 'Study_E4/Nr_1/';
autorun.alg.Sim_ID{2}  = 'Study_E4/Nr_2/';
autorun.alg.Sim_ID{3}  = 'Study_E4/Nr_3/';
autorun.alg.Sim_ID{4}  = 'Study_E4/Nr_4/';
autorun.alg.Sim_ID{5}  = 'Study_E4/Nr_5/';
autorun.alg.Sim_ID{6}  = 'Study_E4/Nr_6/';
autorun.alg.Sim_ID{7}  = 'Study_E4/Nr_7/';
autorun.alg.Sim_ID{8}  = 'Study_E4/Nr_8/';
autorun.alg.Sim_ID{9}  = 'Study_E4/Nr_9/';
autorun.alg.Sim_ID{10} = 'Study_E4/Nr_10/';
autorun.alg.Sim_ID{11} = 'Study_E4/Nr_11/';
autorun.alg.Sim_ID{12} = 'Study_E4/Nr_12/';
autorun.alg.Sim_ID{13} = 'Study_E4/Nr_13/';
autorun.alg.Sim_ID{14} = 'Study_E4/Nr_14/';
autorun.alg.Sim_ID{15} = 'Study_E4/Nr_15/';
autorun.alg.Sim_ID{16} = 'Study_E4/Nr_16/';
autorun.alg.Sim_ID{17} = 'Study_E4/Nr_17/';
autorun.alg.Sim_ID{18} = 'Study_E4/Nr_18/';
autorun.alg.Sim_ID{19} = 'Study_E4/Nr_19/';
autorun.alg.Sim_ID{20} = 'Study_E4/Nr_20/';

% [m/s] velocity of upper body in x1-direction, the disk
autorun.opc.u_up_ini(1:10)            = 2e-2;           
autorun.opc.u_up_ini(11:20)           = 1e-0; 

% [-]       number of discretization cells
autorun.geo.Nx1(1)  = 40+1;
autorun.geo.Nx1(2)  = 80+1;
autorun.geo.Nx1(3)  = 160+1;
autorun.geo.Nx1(4)  = 320+1;
autorun.geo.Nx1(5)  = 640+1;
autorun.geo.Nx1(6)  = 40+1;
autorun.geo.Nx1(7)  = 80+1;
autorun.geo.Nx1(8)  = 160+1;
autorun.geo.Nx1(9)  = 320+1;
autorun.geo.Nx1(10) = 640+1;
autorun.geo.Nx1(11) = 40+1;
autorun.geo.Nx1(12) = 80+1;
autorun.geo.Nx1(13) = 160+1;
autorun.geo.Nx1(14) = 320+1;
autorun.geo.Nx1(15) = 640+1;
autorun.geo.Nx1(16) = 40+1;
autorun.geo.Nx1(17) = 80+1;
autorun.geo.Nx1(18) = 160+1;
autorun.geo.Nx1(19) = 320+1;
autorun.geo.Nx1(20) = 640+1;

% [-]   flag, which discretization scheme is used for the Couette term 
flag_Couette_Discrete_1 = 1; % 1.Order Upwind
flag_Couette_Discrete_2 = 2; % 2.Order QUICK
autorun.alg.FBNS.flag_Couette_Discrete(1:5)   = flag_Couette_Discrete_1;  
autorun.alg.FBNS.flag_Couette_Discrete(6:10)  = flag_Couette_Discrete_2;     
autorun.alg.FBNS.flag_Couette_Discrete(11:15) = flag_Couette_Discrete_1;  
autorun.alg.FBNS.flag_Couette_Discrete(16:20) = flag_Couette_Discrete_2; 
    
autorun.sim.it = 1;
while autorun.sim.it <= autorun.sim.N
    run EHL_01_setup_Study_E4_gen.m
    run EHL_02_mainprocess.m
    autorun.sim.it = autorun.sim.it + 1;
end
