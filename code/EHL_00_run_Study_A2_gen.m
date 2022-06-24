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
autorun.sim.N = 12;      
% Set parameters which are different for the simulations within this study:
% Study and Simualtion Identifiers:
autorun.alg.Sim_ID{1}   = 'Study_A2/Nr_1/';
autorun.alg.Sim_ID{2}   = 'Study_A2/Nr_2/';
autorun.alg.Sim_ID{3}   = 'Study_A2/Nr_3/';
autorun.alg.Sim_ID{4}   = 'Study_A2/Nr_4/';
autorun.alg.Sim_ID{5}   = 'Study_A2/Nr_5/';
autorun.alg.Sim_ID{6}   = 'Study_A2/Nr_6/';
autorun.alg.Sim_ID{7}   = 'Study_A2/Nr_7/';
autorun.alg.Sim_ID{8}   = 'Study_A2/Nr_8/';
autorun.alg.Sim_ID{9}   = 'Study_A2/Nr_9/';
autorun.alg.Sim_ID{10}  = 'Study_A2/Nr_10/';
autorun.alg.Sim_ID{11}  = 'Study_A2/Nr_11/';
autorun.alg.Sim_ID{12}  = 'Study_A2/Nr_12/';

% [m/s] mean velocities for each operating condition, must be larger than 0
u_m = 9.0e-2;
autorun.opc.u_m_ini(1:12)    =  u_m;                                                                                    

% [-] slide-to-roll ratio // -1 for u_low=0 // 0 for u_low=u_up 
srr_1 =  0.0e-1;
srr_2 = -5.0e-1;
autorun.opc.srr_ini(1:12)    = srr_1;       

% [m] microcavity radius (different for fig 5 & 8 in Mourier (2006))
r_1 = 15.5e-6; 
r_2 = 21.5e-6;  
autorun.geo_Ti.r(1:12)       =  r_1;   

% [m] microcavity depth from Mourier (different for simulation & experiment)
d_1 = 7e-6;
d_2 = 175e-9;
d_3 = 1.3e-6;
d_4 = 160e-9;
autorun.geo_Ti.d(1:6)        =  d_1;   
autorun.geo_Ti.d(7:12)       =  d_2;  

% [-] number of discretization points in the x1-direction
Nx1_1 = 33;
Nx1_2 = 65;
Nx1_3 = 129;
autorun.geo.Nx1(1)  = Nx1_1;
autorun.geo.Nx1(2)  = Nx1_2;
autorun.geo.Nx1(3)  = Nx1_3;
autorun.geo.Nx1(4)  = Nx1_1;
autorun.geo.Nx1(5)  = Nx1_2;
autorun.geo.Nx1(6)  = Nx1_3;
autorun.geo.Nx1(7)  = Nx1_1;
autorun.geo.Nx1(8)  = Nx1_2;
autorun.geo.Nx1(9)  = Nx1_3;
autorun.geo.Nx1(10) = Nx1_1;
autorun.geo.Nx1(11) = Nx1_2;
autorun.geo.Nx1(12) = Nx1_3;


% [-]   flag, which discretization scheme is used for the Couette term 
flag_Couette_Discrete_1 = 1; % 1.Order Upwind
flag_Couette_Discrete_2 = 2; % 2.Order QUICK
autorun.alg.FBNS.flag_Couette_Discrete(1:3)     = flag_Couette_Discrete_1;  
autorun.alg.FBNS.flag_Couette_Discrete(4:6)     = flag_Couette_Discrete_2;  
autorun.alg.FBNS.flag_Couette_Discrete(7:9)     = flag_Couette_Discrete_1;  
autorun.alg.FBNS.flag_Couette_Discrete(10:12)   = flag_Couette_Discrete_2; 
    
autorun.sim.it = 1;
while autorun.sim.it <= autorun.sim.N
    run EHL_01_setup_Study_A2_gen.m
    run EHL_02_mainprocess.m
    autorun.sim.it = autorun.sim.it + 1;
end