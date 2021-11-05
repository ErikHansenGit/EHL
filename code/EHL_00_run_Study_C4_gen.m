close all; clc; clearvars;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EHL study A setup and run code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main developer:
% Erik Hansen: erik.hansen@kit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information:
% This script automatically executes all setups and subsequent simulations
% of study C to replicate the setup of the Multiple dimples of
% Woloszynski et al., 2015:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [-] amount of different simulations to be run:
autorun.sim.N = 32;      
% Set parameters which are different for the simulations within this study:
% Study and Simualtion Identifiers:
autorun.alg.Sim_ID{1}  = 'Study_C4/Nr_1/';
autorun.alg.Sim_ID{2}  = 'Study_C4/Nr_2/';
autorun.alg.Sim_ID{3}  = 'Study_C4/Nr_3/';
autorun.alg.Sim_ID{4}  = 'Study_C4/Nr_4/';
autorun.alg.Sim_ID{5}  = 'Study_C4/Nr_5/';
autorun.alg.Sim_ID{6}  = 'Study_C4/Nr_6/';
autorun.alg.Sim_ID{7}  = 'Study_C4/Nr_7/';
autorun.alg.Sim_ID{8}  = 'Study_C4/Nr_8/';

autorun.alg.Sim_ID{9}  = 'Study_C4/Nr_9/';
autorun.alg.Sim_ID{10} = 'Study_C4/Nr_10/';
autorun.alg.Sim_ID{11} = 'Study_C4/Nr_11/';
autorun.alg.Sim_ID{12} = 'Study_C4/Nr_12/';
autorun.alg.Sim_ID{13} = 'Study_C4/Nr_13/';
autorun.alg.Sim_ID{14} = 'Study_C4/Nr_14/';
autorun.alg.Sim_ID{15} = 'Study_C4/Nr_15/';
autorun.alg.Sim_ID{16} = 'Study_C4/Nr_16/';

autorun.alg.Sim_ID{17} = 'Study_C4/Nr_17/';
autorun.alg.Sim_ID{18} = 'Study_C4/Nr_18/';
autorun.alg.Sim_ID{19} = 'Study_C4/Nr_19/';
autorun.alg.Sim_ID{20} = 'Study_C4/Nr_20/';
autorun.alg.Sim_ID{21} = 'Study_C4/Nr_21/';
autorun.alg.Sim_ID{22} = 'Study_C4/Nr_22/';
autorun.alg.Sim_ID{23} = 'Study_C4/Nr_23/';
autorun.alg.Sim_ID{24} = 'Study_C4/Nr_24/';

autorun.alg.Sim_ID{25} = 'Study_C4/Nr_25/';
autorun.alg.Sim_ID{26} = 'Study_C4/Nr_26/';
autorun.alg.Sim_ID{27} = 'Study_C4/Nr_27/';
autorun.alg.Sim_ID{28} = 'Study_C4/Nr_28/';
autorun.alg.Sim_ID{29} = 'Study_C4/Nr_29/';
autorun.alg.Sim_ID{30} = 'Study_C4/Nr_30/';
autorun.alg.Sim_ID{31} = 'Study_C4/Nr_31/';
autorun.alg.Sim_ID{32} = 'Study_C4/Nr_32/';

% [-]   flag whether solid body is elastic on macroscopic scale (true) or rigid (false)
autorun.sld.flag_elastic(1:8)   = false;     
autorun.sld.flag_elastic(9:16)  = true; 
autorun.sld.flag_elastic(17:24) = false; 
autorun.sld.flag_elastic(25:32) = true; 

% [-]   relaxation factor for update of hydrodynamic pressure
autorun.FBNS.alpha_p_nd(1:8)    = 1e-0;         
autorun.FBNS.alpha_p_nd(9:16)   = 5e-1;         
autorun.FBNS.alpha_p_nd(17:24)  = 1e-0;         
autorun.FBNS.alpha_p_nd(25:32)  = 5e-1; 
% [-]   relaxation factor for update of cavity fraction
autorun.FBNS.alpha_thet(1:8)     = 1e-0;         
autorun.FBNS.alpha_thet(9:16)    = 5e-1;         
autorun.FBNS.alpha_thet(17:24)   = 1e-0;         
autorun.FBNS.alpha_thet(25:32)   = 5e-1;        


% [-]       number of dimples
autorun.geo.K_x1(1)            = 1e0;           
autorun.geo.K_x1(2)            = 2e0;          
autorun.geo.K_x1(3)            = 4e0;           
autorun.geo.K_x1(4)            = 8e0;         
autorun.geo.K_x1(5)            = 1e1;          
autorun.geo.K_x1(6)            = 2e1;          
autorun.geo.K_x1(7)            = 4e1;          
autorun.geo.K_x1(8)            = 8e1;           

autorun.geo.K_x1(9)            = 1e0;          
autorun.geo.K_x1(10)           = 2e0;         
autorun.geo.K_x1(11)           = 4e0;           
autorun.geo.K_x1(12)           = 8e0;          
autorun.geo.K_x1(13)           = 1e1;          
autorun.geo.K_x1(14)           = 2e1;         
autorun.geo.K_x1(15)           = 4e1;          
autorun.geo.K_x1(16)           = 8e1;    

autorun.geo.K_x1(17)           = 1e0;           
autorun.geo.K_x1(18)           = 2e0;          
autorun.geo.K_x1(19)           = 4e0;           
autorun.geo.K_x1(20)           = 8e0;         
autorun.geo.K_x1(21)           = 1e1;          
autorun.geo.K_x1(22)           = 2e1;          
autorun.geo.K_x1(23)           = 4e1;          
autorun.geo.K_x1(24)           = 8e1;           

autorun.geo.K_x1(25)           = 1e0;          
autorun.geo.K_x1(26)           = 2e0;         
autorun.geo.K_x1(27)           = 4e0;           
autorun.geo.K_x1(28)           = 8e0;          
autorun.geo.K_x1(29)           = 1e1;          
autorun.geo.K_x1(30)           = 2e1;         
autorun.geo.K_x1(31)           = 4e1;          
autorun.geo.K_x1(32)           = 8e1;  

% [-]   flag, which discretization scheme is used for the Couette term 
flag_Couette_Discrete_1 = 1; % 1.Order Upwind
flag_Couette_Discrete_2 = 2; % 2.Order QUICK
autorun.alg.FBNS.flag_Couette_Discrete(1:16)  = flag_Couette_Discrete_1;  
autorun.alg.FBNS.flag_Couette_Discrete(17:32) = flag_Couette_Discrete_2;     

    
autorun.sim.it = 1;
while autorun.sim.it <= autorun.sim.N
    run EHL_01_setup_Study_C4_gen.m
    run EHL_02_mainprocess.m
    autorun.sim.it = autorun.sim.it + 1;
end