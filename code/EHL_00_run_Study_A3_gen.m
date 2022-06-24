close all; clc; clearvars;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EHL study A setup and run code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main developer:
% Erik Hansen: erik.hansen@kit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information:
% This script automatically executes all setups and subsequent simulations
% of study A to replicate the setup of the Ball on Disc Tribometer of
% Mourier et al., 2006.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [-] amount of different simulations to be run:
autorun.sim.N = 8;      
% Set parameters which are different for the simulations within this study:
% Study and Simualtion Identifiers:
autorun.alg.Sim_ID{1} = 'Study_A3/Nr_1/';
autorun.alg.Sim_ID{2} = 'Study_A3/Nr_2/';
autorun.alg.Sim_ID{3} = 'Study_A3/Nr_3/';
autorun.alg.Sim_ID{4} = 'Study_A3/Nr_4/';
autorun.alg.Sim_ID{5} = 'Study_A3/Nr_5/';
autorun.alg.Sim_ID{6} = 'Study_A3/Nr_6/';
autorun.alg.Sim_ID{7} = 'Study_A3/Nr_7/';
autorun.alg.Sim_ID{8} = 'Study_A3/Nr_8/';

% [m/s] mean velocities for each operating condition, must be larger than 0
u_m = 9.0e-2;
autorun.opc.u_m_ini(1:8)    =  u_m;                                                                                    

% [-] slide-to-roll ratio // -1 for u_low=0 // 0 for u_low=u_up 
srr_1 =  0.0e-1;
srr_2 = -5.0e-1;
autorun.opc.srr_ini(1:2)    = srr_1;                                             
autorun.opc.srr_ini(3:4)    = srr_2;        
autorun.opc.srr_ini(5:6)    = srr_1;           
autorun.opc.srr_ini(7:8)    = srr_2;         

% [m] microcavity radius (different for fig 5 & 8 in Mourier (2006))
r_1 = 15.5e-6; 
r_2 = 21.5e-6;  
autorun.geo_Ti.r(1:2)       =  r_1;   
autorun.geo_Ti.r(3:4)       =  r_2;   
autorun.geo_Ti.r(5:6)       =  r_1;   
autorun.geo_Ti.r(7:8)       =  r_2;   

% [m] microcavity depth from Mourier (different for simulation & experiment)
d_1 = 7e-6;
d_2 = 175e-9;
d_3 = 1.3e-6;
d_4 = 160e-9;
autorun.geo_Ti.d(1)       =  d_1;   
autorun.geo_Ti.d(2)       =  d_2;   
autorun.geo_Ti.d(3)       =  d_3;   
autorun.geo_Ti.d(4)       =  d_4;   
autorun.geo_Ti.d(5)       =  d_1;   
autorun.geo_Ti.d(6)       =  d_2;   
autorun.geo_Ti.d(7)       =  d_3;   
autorun.geo_Ti.d(8)       =  d_4;   

% [-]   flag, which discretization scheme is used for the Couette term 
flag_Couette_Discrete_1 = 1; % 1.Order Upwind
flag_Couette_Discrete_2 = 2; % 2.Order QUICK
autorun.alg.FBNS.flag_Couette_Discrete(1:4) = flag_Couette_Discrete_1;  
autorun.alg.FBNS.flag_Couette_Discrete(5:8) = flag_Couette_Discrete_2;     
    
autorun.sim.it = 1;
while autorun.sim.it <= autorun.sim.N
    run EHL_01_setup_Study_A3_gen.m
    run EHL_02_mainprocess.m
    autorun.sim.it = autorun.sim.it + 1;
end