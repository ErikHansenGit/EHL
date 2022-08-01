close all; clc; clearvars -except autorun;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EHL mainprocess code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main developer:
% Erik Hansen: erik.hansen@kit.edu
% Further developers who contributed to the code:
% Altay Kacan, 08.01.2021: Implementation of the unsteady terms in the
% Reynolds equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information:
% This code solves the (Unsteady) Reynolds equation under consideration of a
% mass conserving cavitation model
% Furthermore, compressibility, piezoviscous and shear thinning (to
% a certain degree) effects and
% elastic surface deformation due to the fluid pressure can be modeled. 
% Either a constant imposed rigid body displacement or a constant imposed
% normal load can be specified.
%
% The derivation of the Reynolds equation with mass-conserving cavitation
% through the cavity fraction is described by:
% Giacopini, M., Fowell, M. T., Dini, D., & Strozzi, A. (2010). 
% A mass-conserving complementarity formulation to study lubricant films in the presence of cavitation. 
% Journal of tribology, 132(4).
% 
% The computation of the hydrodynamic pressure and cavity fraction distribution for
% a given geometry with the Fischer-Burmeister-Newton-Schur (FBSN) algorithm
% is done according to:
% Woloszynski, T., Podsiadlo, P., & Stachowiak, G. W. (2015). 
% Efficient solution to the cavitation problem in hydrodynamic lubrication.
% Tribology Letters, 58(1), 18.
% The this code was largely inspired by the FBNS implementation of:
% Codrignani, A. R. (2019). Numerical representation of a pin-on-disc tribometer for the investigation of textured surfaces (Vol. 6). KIT Scientific Publishing.
%
% Even though this solver does not use multigrid methods, the incorporation
% of the elastic deformation into the FBNS algorithm by reducing the Kernel
% influence within the Jacobian of the FBNS algorithm to the five main
% diagonals and using the Roelands and Dowson-Higginson relations for the 
% dependence of dynamic viscosity and density of the liquid was largely inspired by:
% Venner, C. H., & Lubrecht, A. A. (Eds.). (2000). Multi-level methods in lubrication. Elsevier.
%
% The consideration of unsteady phenomena is based on the preliminary work of 
% Altay Kacan. For details, see his Bachelor's Thesis "Unsteady Numerical
% Simulation of the Pressure Distribution in a Ball-on-Disc Tribometer",
% 2021, Karlsruhe Institute of Technology
%
% The unsteady Reynolds equation is discretized with a 2D Version of the 
% Finite Volume Method (FVM).
% The Poiseuille term is discretized with a second order
% central differential scheme. The Couette Term can be discretized with
% different first, second or thrid order upwind schemes. The time dependent 
% is discretized with the first Order Euler Implicit scheme. 
% Integrals over the cell or cell boundaries are
% approximated through the second order midpoint rule, thus limiting the
% maximum reachable discretization order to second order.
%
% The elastic deformation is considered through the elastic half-space
% model when distinctive constant pressures over single rectangles are assumed
% and the Boundary Lement Method (BEM) is applied.
% The Kernel is constructed as explained by 
% Johnson, K. L., 2004. Contact mechanics. Cambridge: Cambridge University Press
% in equation (3.25) on P.54
% The elastic deformation is computed with a linear convolution of the
% Kernel with the hydrodynamic pressure in Fourier space.
%
% The initial guess for the pressure disttribution can either be taken to
% be ambient pressure or as the solution of the dry contact problem. The
% code which computes the sokution of the dry contact problem is taken from
% https://github.com/ErikHansenGit/Contact_elastic_half-space
%
% The adjustment of the rigid body displacement to meet an imposed normal
% force is realized with a PID-controller as motivated and explained by:
% Wang, Yuechang, Ying Liu and Yuming Wang, 2019. A method for improving the 
% capability of convergence of numerical lubrication simulation by using the 
% PID controller. In: IFToMM World Congress on Mechanism and Machine Science. 
% Springer. 2019. p. 3845–3854. IFToMM World Congress on Mechanism and Machine 
% Science
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read in setup data:
% Input path:
input_path  =              './../data/EHL_02_mainprocess/Input/';
% Load everything:
load(fullfile(input_path,'alg.mat'));
load(fullfile(input_path,'fld.mat'));
load(fullfile(input_path,'sld.mat'));
load(fullfile(input_path,'geo.mat'));
load(fullfile(input_path,'opc.mat'));
if alg.flag_unsteady == 1 
    load(fullfile(input_path,'geo_Ti.mat'));
end
% Create Output directories
output_main_path = char(compose('./../data/%s/Output/',alg.Sim_ID));
warning('off','MATLAB:MKDIR:DirectoryExists')
mkdir (output_main_path)
% Create results folder:
output_result_path = fullfile(output_main_path,'Result/');
mkdir (output_result_path)
rmdir(output_result_path, 's')        % remove old results if they exist
mkdir (output_result_path)
% Save used input 
output_input_path = fullfile(output_main_path,'Used_input/');
mkdir (output_input_path)
rmdir(output_input_path, 's')         % remove old inputs if they exist
mkdir (output_input_path)

% Save used input: 
save (fullfile(output_input_path,'fld.mat'),'fld');
save (fullfile(output_input_path,'sld.mat'),'sld');
save (fullfile(output_input_path,'geo.mat'),'geo');
save (fullfile(output_input_path,'opc.mat'),'opc');
save (fullfile(output_input_path,'alg.mat'),'alg');

if alg.flag_unsteady
    
    save (fullfile(output_input_path,'geo_Ti.mat'),'geo_Ti');
    clear geo_Ti;
    
    for it_OC = 1:opc.N
        for it_time = 1:opc.N_t(it_OC)
            % Define path from where to load and save the time dependent geometry
            % variation
            input_path_OC   = fullfile(input_path,sprintf('h_time/OC_%i/',it_OC));        
            input_path_time = fullfile(input_path_OC,sprintf('Time_%i/',it_time));
            
            % Load data:
            load(fullfile(input_path_time,'slc.mat'));

            % Save used input:
            % Define paths:  
            output_input_path_OC   = fullfile(output_input_path,sprintf('h_time/OC_%i/',it_OC));
            output_input_path_time = fullfile(output_input_path_OC,sprintf('Time_%i/',it_time));
            % Create paths:            
            mkdir (output_input_path_OC)
            mkdir (output_input_path_time)
            % Save read in data:
            save(fullfile(output_input_path_time,'slc.mat'),'slc');
        end
    end
    clear slc;
end


%% Start actual simulation:
% Apply performance settings:
% Limit maximum number of computational threads if desired:
if alg.flag_limit_num_comp_threads
    maxNumCompThreads(alg.max_num_comp_threads)
end
% Activate profile viewer:
if alg.flag_profile_viewer
    profile on
end

fprintf('\n--------------------------------');
fprintf('\n<strong>Simulation started!</strong>');
fprintf('\nSimulation ID: ''%s''',alg.Sim_ID);
fprintf('\n--------------------------------\n');

% Start timer:
time_start = tic;

%% Get remaining information and initialize:
h.Nx1 = geo.Nx1;
h.Nx2 = geo.Nx2;
h.N = geo.N;
h.x1 = geo.x1;
h.x2 = geo.x2;
h.dx1 = geo.dx1;
h.dx2 = geo.dx2;
h.h_g_ma = geo.h_g_ma;
clear geo;

if opc.flag_imposed == 1
    h.h_d_ma = opc.h_d*ones(alg.it_tot_max,1);
elseif opc.flag_imposed == 2
    h.h_d_ma = zeros(alg.it_tot_max,1);
    h.h_d_ma(1) = alg.load.h_d_ma_ini;    
end

h.h_el_hd_ma    = zeros(h.Nx1,h.Nx2);
h.h_ma          = zeros(h.Nx1,h.Nx2);
h.h_m           = zeros(h.Nx1,h.Nx2);

% Compute linear indices
% MATLAB first linearly indices along columns -> i -> x1 then along rows -> j -> x2:
[i_matr, j_matr]    = ndgrid(1:h.Nx1,1:h.Nx2);
h.i_lin             = zeros(h.Nx1,h.Nx2);
h.i_lin(:,:)        = sub2ind([h.Nx1 h.Nx2], i_matr(1:h.Nx1,1:h.Nx2), j_matr(1:h.Nx1,1:h.Nx2)); % [-]   linear indizes of the discretization points

% Kernel
if sld.flag_elastic
    % Compute kernel of elastic deformations:
    [h.fft2_Kernel,h.Kernel]  = construct_linear_Kernel(h.Nx1,h.Nx2,h.dx1,h.dx2,sld.E_dash);
end

%% Initialisations:
% Initial guess of the hydrodynamic pressure field:
if alg.FBNS.flag_p_hd_ini==1    
    alg.FBNS.p_hd_ini.p_hd   = fld.p_amb*ones(h.Nx1,h.Nx2);      % [Pa]  initial guess for the hydrodynamic pressure

elseif alg.FBNS.flag_p_hd_ini==2
    p_hd_ini.z               = -h.h_g_ma;                        % [-]   geometry profile as an input for the contact algorithm
    p_hd_ini.h_s             = h.h_d_ma(1) - fld.h_min;          % [-]   seperating distance between upper and lower reference planes as an input for the contact algorithm
    p_hd_ini.H               = Inf;                              % [Pa]  upper limit for the initial guess of the contact pressure
    p_hd_ini.p_con_ini       = zeros(h.Nx1,h.Nx2);               % [Pa]  initial guess for contact algorithm
    [alg.FBNS.p_hd_ini.p_hd,~,res.FBNS.p_hd_ini]= elpl_contact_pressure_akchurin_linear(fld.p_amb,p_hd_ini.H,...
        p_hd_ini.p_con_ini,p_hd_ini.z,p_hd_ini.h_s,h.Nx1,h.Nx2,h.fft2_Kernel,...
        alg.FBNS.p_hd_ini.toll,alg.FBNS.p_hd_ini.it_max,alg.FBNS.ref.h);
    alg.FBNS.p_hd_ini.it      = length(res.FBNS.p_hd_ini);       % [-]   number of contact algorithm iterations
        
    fprintf('\n--------------------------------');
    fprintf('\nInitial pressure was computed');
    fprintf('\nthrough dry contact problem');
    fprintf('\nit: con         = %i',alg.FBNS.p_hd_ini.it);
    fprintf('\nres: con        = %e',res.FBNS.p_hd_ini(alg.FBNS.p_hd_ini.it));
    fprintf('\n--------------------------------\n');
end

% State that no previous results exist for the load balance PID controller:
PID.flag_ini = false;

%% Wrappers
% Operating condition wrapper
for it_OC= 1:opc.N                                                      % [-]       operating condition iteration counter
    
    alg.time.N = opc.N_t(it_OC);                                        % [-]   amount of discrete time steps 
    if alg.flag_unsteady
        alg.time.delta_t = opc.delta_t(it_OC);                          % [s]  time step size
    end
    
    if it_OC > 1
        if opc.flag_imposed == 1
            h.h_d_ma = opc.h_d*ones(alg.it_tot_max,1);
        elseif opc.flag_imposed == 2
            h.h_d_ma = zeros(alg.it_tot_max,1);
            h.h_d_ma(1) = alg.load.h_d_ma_ini;    
        end

        h.h_el_hd_ma    = zeros(h.Nx1,h.Nx2);
        h.h_ma          = zeros(h.Nx1,h.Nx2);
        h.h_m           = zeros(h.Nx1,h.Nx2);
    end
    
    % Time iteration wrapper
    for it_time = 1:alg.time.N                                          % [-]       time iteration counter

        %% Setting initial guesses:
        if it_time == 1 
            sol.p_hd                = alg.FBNS.p_hd_ini.p_hd;           % [Pa]      solution of the hydrodynamic pressure field
            sol.thet                = zeros(h.Nx1,h.Nx2);               % [-]       solution of the cavity fraction
        else
            % Use solution of previous time step as initial guess
            sol.p_hd                = prev.sol.p_hd;                    % [Pa]      solution of the hydrodynamic pressure field
            sol.thet                = prev.sol.thet;                    % [-]       solution of the cavity fraction
        end

        %% Initialisations:
        alg.it_tot              = 0;                                % [-]       total iteration counter
        alg.time.it             = it_time;                          % [-]       time iteration counter
        alg.opc.it              = it_OC;                            % [-]       time iteration counter

        if opc.flag_imposed == 2
            alg.load.F_N_aim    = opc.F_N;                          % [N]       desired load force
        end
        
        if it_time > 1
            % Reuse previous rigid body displacement:
            if opc.flag_imposed == 1
                h.h_d_ma    = opc.h_d*ones(alg.it_tot_max,1);
            elseif opc.flag_imposed == 2
                h.h_d_ma    = zeros(alg.it_tot_max,1);
                h.h_d_ma(1) = prev.h.h_d_ma;    
            end
            % Reset:
            h.h_el_hd_ma    = zeros(h.Nx1,h.Nx2);
            h.h_ma          = zeros(h.Nx1,h.Nx2);
            h.h_m           = zeros(h.Nx1,h.Nx2);
        end
        
        % Results of previous time steps (only really used if alg.time.it > 1)
        if alg.time.it == 1 
            prev.prop_nd.rho = zeros(h.Nx1,h.Nx2);
            prev.h_nd.h_m    = zeros(h.Nx1,h.Nx2);
            prev.sol.p_hd    = zeros(h.Nx1,h.Nx2);
            prev.sol.thet    = zeros(h.Nx1,h.Nx2);    
            prev.G           = zeros(h.N,1);
            prev.F           = zeros(h.N,1);         
            prev.h.h_d_ma    = 0;          
        end
        
        % Clock time evaluations:
        alg.FBNS.eval_time.tot       = zeros(1,1);                      % [s]       time to perform the FBNS algorithm
        alg.FBNS.eval_time.it        = zeros(alg.it_tot_max,1);         % [s]       time to perform one iteration in the FBNS algorithm
        alg.FBNS.eval_time.solve_LES = zeros(alg.it_tot_max,1);         % [s]       time to solve the linear equation system in the FBNS algorithm


        % Residuals:
        res.FBNS.FBNS               = zeros(alg.it_tot_max,1);         % [-]       resulting residual of the FBNS algorithm
        res.FBNS.delta_p_nd_mean    = zeros(alg.it_tot_max,1);         % [-]
        res.FBNS.delta_thet_mean    = zeros(alg.it_tot_max,1);         % [-]
        res.FBNS.delta_p_nd_max     = zeros(alg.it_tot_max,1);         % [-]
        res.FBNS.delta_thet_max     = zeros(alg.it_tot_max,1);         % [-]
        res.FBNS.delta_G_max        = zeros(alg.it_tot_max,1);         % [-]
        res.FBNS.delta_F_max        = zeros(alg.it_tot_max,1);         % [-]
        res.FBNS.G_max              = zeros(alg.it_tot_max,1);         % [-]
        res.FBNS.F_max              = zeros(alg.it_tot_max,1);         % [-]
        res.load.F_N                = zeros(alg.it_tot_max,1);         % [-

        %% Dimensional quantities:
        prop.mu_l   = fld.mu_0*ones(h.Nx1,h.Nx2);           % [Pas]     dynamic viscosity of liquid lubricant
        prop.rho_l  = fld.rho_0*ones(h.Nx1,h.Nx2);          % [kg/m^3]  density of liquid lubricant
        prop.u_m    = opc.u_m(it_OC);                       % [m/s]     mean velocity
        prop.u_r    = opc.u_r(it_OC);                       % [m/s]     relative velocity

        if alg.flag_unsteady
            % Define path from where to load the time dependent geometry
            % variation
            % Define paths:  
            output_input_path_OC   = fullfile(output_input_path,sprintf('h_time/OC_%i/',it_OC));
            output_input_path_time = fullfile(output_input_path_OC,sprintf('Time_%i/',it_time));
           
            % Load data:
            load(fullfile(output_input_path_time,'slc.mat'));
            
            % Apply geometry variation:
            h.h_t_ma    = slc.h_t_ma;           % [m]       time dependent gap height variation
            alg.time.t  = slc.t;                % [s]       corresponding time
            
        end

        %% Reference quantities:
        ref.x1  = alg.FBNS.ref.x1;
        ref.x2  = alg.FBNS.ref.x2;
        ref.h   = alg.FBNS.ref.h;
        ref.mu  = alg.FBNS.ref.mu;
        ref.rho = alg.FBNS.ref.rho;
        ref.u   = alg.FBNS.ref.u(it_OC);
        ref.p   = alg.FBNS.ref.p;
        if alg.flag_unsteady
            ref.t   = alg.FBNS.ref.t(it_OC);
        end
        if opc.flag_imposed == 2
            ref.F_N   = alg.load.ref.F_N;
        end

        %% Non-dimensional quantities:
        h_nd.Nx1    = h.Nx1;
        h_nd.Nx2    = h.Nx2;
        h_nd.N      = h.N;
        h_nd.i_lin  = h.i_lin;
        h_nd.dx1    = h.dx1             /ref.x1;
        h_nd.dx2    = h.dx2             /ref.x2;
        h_nd.h_d_ma = h.h_d_ma          /ref.h;
        h_nd.h_m    = h.h_m             /ref.h;
        h_nd.x1     = h.x1              /ref.x1;
        h_nd.x2     = h.x2              /ref.x2;
        if sld.flag_elastic  
            h_nd.Kernel = h.Kernel*ref.p/ref.h;
        end
        if alg.flag_unsteady
            h_nd.h_t_ma = h.h_t_ma/ref.h;   

            alg.time.t_nd      = alg.time.t/ref.t;
            alg.time.delta_t_nd = alg.time.delta_t/ref.t;                    
        end

        prop_nd.rho = prop.rho_l    /ref.rho;
        prop_nd.mu  = prop.mu_l     /ref.mu;
        prop_nd.u   = prop.u_m      /ref.u;

        fld_nd.flag_p_bound = fld.flag_p_bound;
        fld_nd.p_amb    = fld.p_amb   /ref.p;
        fld_nd.p_cav    = fld.p_cav   /ref.p;
        fld_nd.p_S_side = fld.p_S_side/ref.p;
        fld_nd.p_W_side = fld.p_W_side/ref.p;
        fld_nd.p_E_side = fld.p_E_side/ref.p;
        fld_nd.p_N_side = fld.p_N_side/ref.p;

        fld_nd.flag_thet_bound = fld.flag_thet_bound;
        fld_nd.thet_S_side = fld.thet_S_side;
        fld_nd.thet_W_side = fld.thet_W_side;
        fld_nd.thet_E_side = fld.thet_E_side;
        fld_nd.thet_N_side = fld.thet_N_side;

        %% Solve equation system:
        time_FBNS = tic;
        [prop,prop_nd,h_nd,alg,res,sol,h,prev,PID] = FBNS_nd(prop,prop_nd,h_nd,fld,fld_nd,ref,alg,res,sol,sld,h,opc,time_start,prev,PID);
        alg.FBNS.eval_time.tot = toc(time_FBNS);

        % Compute shear stresses on upper and lower surface:
        [sol.tau_hd_low,sol.tau_hd_up,~] = compute_tau_hd(h,sol,prop); 

        %% Stribeck data for each time point:
        sol.F_T_up        = -sum( sol.tau_hd_up( :)*h.dx1*h.dx2);       % [N]   friction force in positive x1-direction acting from the fluid upon the upper surface 
        sol.F_T_low       = -sum(-sol.tau_hd_low(:)*h.dx1*h.dx2);       % [N]   friction force in positive x1-direction acting from the fluid upon the lower surface 
        sol.C_f_up        = abs(sol.F_T_up /sol.F_N);                   % [-]   friction coefficient if measured on upper surface
        sol.C_f_low       = abs(sol.F_T_low/sol.F_N);                   % [-]   friction coefficient if measured on lower surface
        sol.h_m_min       = min(h.h_m(:));                              % [m]   minimum adjusted macroscopic gap height

        % Resize arrays:
        h.h_d_ma      = h.h_d_ma(1:alg.it_tot,1); 
        h_nd.h_d_ma   = h_nd.h_d_ma(1:alg.it_tot,1);
               
        res.FBNS.FBNS               = res.FBNS.FBNS(1:alg.it_tot,1);
        res.FBNS.delta_p_nd_mean    = res.FBNS.delta_p_nd_mean(1:alg.it_tot,1);
        res.FBNS.delta_thet_mean    = res.FBNS.delta_thet_mean(1:alg.it_tot,1); 
        res.FBNS.delta_p_nd_max     = res.FBNS.delta_p_nd_max(1:alg.it_tot,1);
        res.FBNS.delta_thet_max     = res.FBNS.delta_thet_max(1:alg.it_tot,1); 
        res.FBNS.delta_G_max        = res.FBNS.delta_G_max(1:alg.it_tot,1);
        res.FBNS.delta_F_max        = res.FBNS.delta_F_max(1:alg.it_tot,1); 
        res.FBNS.G_max              = res.FBNS.G_max(1:alg.it_tot,1);
        res.FBNS.F_max              = res.FBNS.F_max(1:alg.it_tot,1); 

        res.load.F_N                = res.load.F_N(1:alg.it_tot,1); 

        alg.FBNS.eval_time.it        = alg.FBNS.eval_time.solve_LES(1:alg.it_tot,1);
        alg.FBNS.eval_time.solve_LES = alg.FBNS.eval_time.solve_LES(1:alg.it_tot,1);

        % Consolidate information:
        alg.FBNS.eval_time.solve_LES_tot = sum(alg.FBNS.eval_time.solve_LES(:));      % [s]   total time to solve the linear equation system in the FBNS algorithm
        alg.simulation_total_runtime     = toc(time_start);                           % [s]   total runtime of simulation
        
        % Save results of each operating condition and time point:       
        sub_result_path = sprintf('OC_%i/',it_OC);
        sub_sub_result_path = sprintf('Time_%i/',alg.time.it);
        output_sub_result_path = fullfile(output_result_path,sub_result_path,sub_sub_result_path);
        mkdir (output_sub_result_path)
        save(fullfile(output_sub_result_path,'sol.mat'),'sol');
        save(fullfile(output_sub_result_path,'prop.mat'),'prop');
        save(fullfile(output_sub_result_path,'h.mat'),'h');
        save(fullfile(output_sub_result_path,'alg.mat'),'alg');
        save(fullfile(output_sub_result_path,'res.mat'),'res');
        save(fullfile(output_sub_result_path,'ref.mat'),'ref');

        fprintf('\n--------------------------------');
        fprintf('\n<strong>FBNS algorithm finished!</strong>');
        fprintf('\nalg.it_tot        = %i',alg.it_tot);
        fprintf('\nFBNS algorithm took %fs',alg.FBNS.eval_time.tot);
        fprintf('\n--------------------------------');

    end
end
        
fprintf('\n--------------------------------');
fprintf('\n<strong>Simulation finished!</strong>');
fprintf('\nSimulation ID: ''%s''',alg.Sim_ID);
fprintf('\n--------------------------------\n');

% Show results:
% Detailed call structure and function performance:
if alg.flag_profile_viewer
    profile viewer
end


function [prop,prop_nd,h_nd,alg,res,sol,h,prev,PID] = FBNS_nd(prop,prop_nd,h_nd,fld,fld_nd,ref,alg,res,sol,sld,h,opc,time_start,prev,PID)
% Prepare for iterations:
p_red_nd    = zeros(h_nd.N,1); 
p_red_nd(:) = (sol.p_hd(:) - fld.p_cav)./ref.p; % [-] hydrodynamic pressure minus cavitation pressure in a vector in order of the linear indizes
thet_red    = zeros(h_nd.N,1);
thet_red(:) = sol.thet(:);                  % [-]  cavity fraction in a vector in order of the linear indizes
if alg.time.it == 1
    G_old   = ones(h_nd.N,1);               % [-]  previous residual of Reynolds equation
    F_old   = ones(h_nd.N,1);               % [-]  previous residual of Fischer-Burmeister equation
else
    G_old   = prev.G;
    F_old   = prev.F;
end
while alg.it_tot == 0 || (      (       res.FBNS.FBNS(alg.it_tot,1)  > alg.FBNS.toll  ...
                                 || abs(res.load.F_N( alg.it_tot,1)) > alg.load.toll) ...
                           &&   alg.it_tot < alg.it_tot_max)
    % Update iteration counters:
    alg.it_tot  = alg.it_tot  + 1;
    
    % Start iteration timer:
    time_it = tic;    
    
    % Compute necessary coefficients:
    if sld.flag_elastic
        % Compute elastic gap deformation on macroscopic scale:    
        [h.h_el_hd_ma] = compute_h_el(sol.p_hd,h.Nx1,h.Nx2,h.fft2_Kernel);   
    end
    % Update gap height on macroscopic scale
    h.h_ma      = h.h_d_ma(alg.it_tot) + h.h_g_ma + h.h_el_hd_ma;
    % Add unsteady geometry variation
    if alg.flag_unsteady
        h.h_ma = h.h_t_ma + h.h_ma;
    end
    % Adjust gap height on macroscopic scale to prevent it from becoming zero:
    % That way, J is better conditioned and the gap height cannot become negative
    h.h_m                   = h.h_ma;
    h.h_m(h.h_m<fld.h_min)  = fld.h_min;
    h_nd.h_m                = h.h_m/ref.h;
    
    % Compute properties:
    prop.rho_l  = compute_rho_l(sol.p_hd,fld,h);                        % [kg/m^3]  density of liquid lubricant
    prop_nd.rho = prop.rho_l    /ref.rho;

    prop        = compute_mu_l(sol,fld,prop,h);                         % [Pas]     dynamic viscosity of liquid lubricant  
    prop_nd.mu  = prop.mu_l     /ref.mu; 

    
    % Construction of A, B, c and J_Gp_nd:
    [A,B,c,J_Gp_nd] = FBNS_construct_G_matr(h_nd,prop_nd,ref,alg,fld_nd,sld,prev);
    
    % Capture J_Gthet:
    J_Gthet = B;

    % Compute residuals and center diagonals of J_Fp_nd and J_Fthet:
    [G,F,J_Fp_nd_C,J_Fthet_C] = FBNS_compute_res_nd(A,B,c,p_red_nd,thet_red,h_nd,fld_nd);

    % Swap columns of the Jacobians:
    [A_G,B_G,A_F,B_F,sw_cond] = FBNS_swapping(J_Gp_nd,J_Gthet,J_Fp_nd_C,J_Fthet_C,h_nd);
   
    % Solve LES:
    time_solve_LES = tic;    
    delta_b     = (B_G - A_G*(A_F\B_F))\(-G + A_G*(A_F\F));
    delta_a     = A_F\(-F - B_F*delta_b);
    alg.FBNS.eval_time.solve_LES(alg.it_tot,1) = toc(time_solve_LES);
    
    % Determine updates for p_nd_red and thet:
    delta_p_nd      = delta_a;
    delta_thet      = delta_b;    
    % Swapping again based on the same condition restores the initial order
    delta_p_nd(sw_cond) = delta_b(sw_cond); 
    delta_thet(sw_cond) = delta_a(sw_cond);
    
    % Determine relaxation coefficients:
    alpha_p_nd = alg.FBNS.alpha_p_nd;
    alpha_thet = alg.FBNS.alpha_thet;

    % Update p_red_nd and thet_red:
    p_red_nd       = p_red_nd     + alpha_p_nd*delta_p_nd;
    thet_red       = thet_red     + alpha_thet*delta_thet;

    % Enforce limits:
    % Truncate updated values below 0:
    if alg.FBNS.flag_enforce_limit_p_nd_red_min
        p_red_nd((p_red_nd<0)) = 0;
    end
    if alg.FBNS.flag_enforce_limit_thet_red_min
        thet_red((thet_red<0)) = 0;
    end
    % Truncate updated values above threshold:
    if alg.FBNS.flag_enforce_limit_p_red_max 
    	p_red_nd((p_red_nd>alg.FBNS.limit_p_red_max/ref.p)) = alg.FBNS.limit_p_red_max/ref.p;
    end
    if alg.FBNS.flag_enforce_limit_thet_red_max
    	thet_red((thet_red>1)) = 1;
    end
    
    % Determine dimensional reduced pressure:
    p_red = p_red_nd*ref.p;
    
    % Compute real hydrodynamic pressure field:
    sol.p_hd(:) = p_red(:) + fld.p_cav; 
    % Determine cavity fraction field:
    sol.thet(:) = thet_red(:);
    
    % Compute resulting normal load force:
    sol.F_N     = sum(sum((sol.p_hd - fld.p_amb)*h.dx1*h.dx2));   % [N]   load force in negative x3-direction
    
    % Compute residuals:
    [G_old,F_old,res] = FBNS_determine_residual(res,alg,delta_p_nd,delta_thet,G,F,G_old,F_old,sol,ref,opc);
    
    % Adjust rigid body displacement such that load balance equation is met
    % if desired:
    if (opc.flag_imposed == 2) && (alg.it_tot < alg.it_tot_max)
        if (res.FBNS.FBNS(alg.it_tot,1) > alg.FBNS.toll) || (abs(res.load.F_N( alg.it_tot,1)) > alg.load.toll)
                       
            % Determine incremental PID-Controller values:
            PID.e_k =  res.load.F_N(alg.it_tot,1);
                       
            if ~PID.flag_ini
                PID.e_km1 = 0;
                PID.e_km2 = 0;
            end
                
            % Apply PID controller to adjust non-dimensional rigid body displacement:
            h_nd.h_d_ma(alg.it_tot+1) = alg.load.K_P*(PID.e_k - PID.e_km1) + h_nd.h_d_ma(alg.it_tot);
            h_nd.h_d_ma(alg.it_tot+1) = alg.load.K_I*PID.e_k + h_nd.h_d_ma(alg.it_tot+1);
            h_nd.h_d_ma(alg.it_tot+1) = alg.load.K_D*(PID.e_k - 2 * PID.e_km1 + PID.e_km2) + h_nd.h_d_ma(alg.it_tot+1);
            % Compute dimensional rigid body displacement:
            h.h_d_ma(alg.it_tot+1)  = h_nd.h_d_ma(alg.it_tot+1)*ref.h;
            
            % Save information for next iteration:
            PID.e_km2 = PID.e_km1;
            PID.e_km1 = PID.e_k;
            PID.flag_ini = true;
        end       
    end
    
    % Determine iteration step time:
    alg.FBNS.eval_time.it(alg.it_tot,1) = toc(time_it);
    
    % Print iteration information:
    if alg.FBNS.flag_print_it == 1
        fprintf('\n--------------------------------');
        fprintf('\nSimulation  ''%s''',alg.Sim_ID);
        fprintf('\nis now running for %fs',toc(time_start));
        fprintf('\n<strong>it: opc         = %i/%i</strong>',alg.opc.it,opc.N);
        fprintf('\nu_m             = %fm/s',opc.u_m(alg.opc.it));
        fprintf('\nu_r             = %fm/s',opc.u_r(alg.opc.it));
        if alg.flag_unsteady
            fprintf('\nUnsteady simulation');
            fprintf('\n<strong>it: time        = %i/%i</strong>\n',alg.time.it,alg.time.N);
        else
            fprintf('\nSteady simulation\n\n');
        end
        fprintf('\n<strong>it: tot         = %i</strong>',alg.it_tot);
        fprintf('\n<strong>res: FBNS       = %e</strong>',res.FBNS.FBNS(alg.it_tot,1));
        fprintf('\nres: d_p_nd     = %e',res.FBNS.delta_p_nd_mean(alg.it_tot,1));
        fprintf('\nres: d_thet     = %e',res.FBNS.delta_thet_mean(alg.it_tot,1));
        fprintf('\nres: d_p_nd_max = %e',res.FBNS.delta_p_nd_max(alg.it_tot,1));
        fprintf('\nres: d_thet_max = %e',res.FBNS.delta_thet_max(alg.it_tot,1));
        fprintf('\nres: d_G_max    = %e',res.FBNS.delta_G_max(alg.it_tot,1));
        fprintf('\nres: d_F_max    = %e',res.FBNS.delta_F_max(alg.it_tot,1));
        fprintf('\nres: G_max      = %e',res.FBNS.G_max(alg.it_tot,1));
        fprintf('\nres: F_max      = %e',res.FBNS.F_max(alg.it_tot,1));
        if opc.flag_imposed == 2
            fprintf('\n<strong>res: F_N        = %e</strong>',res.load.F_N(alg.it_tot,1));
            fprintf('\nh: h_d_ma       = %em',h.h_d_ma(alg.it_tot,1));
        
        end
        fprintf('\n--------------------------------\n');
    end
    
end
% Save results needed for next time step:
if alg.flag_unsteady && (alg.time.it < alg.time.N)
    prev.prop_nd.rho    = prop_nd.rho;     
    prev.h_nd.h_m       = h_nd.h_m;   
    prev.sol.p_hd       = sol.p_hd;
    prev.sol.thet       = sol.thet; 
    
    prev.G              = G_old;
    prev.F              = F_old; 

    prev.h.h_d_ma       = h.h_d_ma(alg.it_tot);  
end
end

function [A,B,c,J_Gp_nd] = FBNS_construct_G_matr(h_nd,prop_nd,ref,alg,fld_nd,sld,prev)
% This function constructs the vectors and matrices needed to compute the
% residual of the dimensionless Reynolds equation G = A*p_nd_red + B*thet_red + c and the
% matrix J_Gp_nd needed to determine the update out of the residual
% The boundary conditions of thet_red are later incorporated in F, J_Fp_nd_C and J_Fthet_C and the boundary
% conditions of p_nd_red are incorporated in G, A, B and c
%
% Construction of sparse A:
% Discretization of Poisseuille terms with CENTRAL scheme
%
% Construction of sparse B:
% Discretization of Couette terms with UPWIND scheme
%
% Construction of c:
% Discretization of Couette terms with UPWIND scheme adn time dependend
% terms with Euler IMPLICIT scheme
%
% Construction of sparse J_Gp_nd:
% J_Gp_nd consists of the five main diagonals of A + A_Co_h + A_Ti_h
% Discretization of Poiseuille terms with CENTRAL scheme and of Couette terms with UPWIND scheme
% Time discretization according to Euler IMPLICIT scheme
% Reason:
% When calculating the Jacobian for the Reynolds Equation, the terms that 
% contain h on its own (without being multiplied by theta or pressure)in the discretized equation
% require the substitution of the pressure dependence of h in order to
% give accurate results. Since the elastic deformation at a given grid
% point depends on the pressures of ALL other grid points a simplification
% is necessary because the Jacobian would not be sparse. That is why only the 
% 5 main diagonals are taken of the otherwise non-sparse Jacobian.
% 
% For the steady case this is the A_Co_h matrix. It contains the coefficients
% of the pressure terms that show up after the subsitution of h.
%
% For the unsteady calculations the Jacobian of the Reynolds Equation has an extra term that
% has h on its own which shows up due to the extra coefficients (A_Ti_h).
%
% Generic 2nd Order Upwind Discretization schemes for the Couette term can
% be considered by adjusting alg.FBNS.flag_Couette_Discrete
%
%
%
% Shorter definitions:
N       = h_nd.N;
Nx1     = h_nd.Nx1;
Nx2     = h_nd.Nx2;
i_lin   = h_nd.i_lin;
dx1_nd  = h_nd.dx1;
dx2_nd  = h_nd.dx2;
x_rat   = ref.x1/ref.x2;

% Inside domain excluding boundaries:
N_in        = (Nx1 - 2)*(Nx2 - 2);      % [-]   number of points within domain
i_in        = 2:Nx1-1;                  % [-]   regular indices in x1 direction
j_in        = 2:Nx2-1;                  % [-]   regular indices in x2 direction
i_lin_in    = i_lin(i_in,j_in);               % [-]   linear indices of the points not on the boundary

% Cells inside domain but without a WestWest neighbour:
N_nWW       = (Nx2 - 2);                % [-]   number of points 
i_nWW       = 2;                        % [-]   regular indices in x1 direction
j_nWW       = 2:Nx2-1;                  % [-]   regular indices in x2 direction
i_lin_nWW   = i_lin(i_nWW,j_nWW);       % [-]   linear indices of the points

% Cells inside domain and with a WestWest neighbour:
N_wWW        = (Nx1 - 3)*(Nx2 - 2);     % [-]   number of points
i_wWW        = 3:Nx1-1;                 % [-]   regular indices in x1 direction
j_wWW        = 2:Nx2-1;                 % [-]   regular indices in x2 direction
i_lin_wWW    = i_lin(i_wWW,j_wWW);      % [-]   linear indices of the points


% Poiseuille term:
% Disretized with central scheme:
% Obtain coefficients at the cell centers:
xi_nd_Po_C      = prop_nd.rho.*(h_nd.h_m   .^3)./prop_nd.mu;

% Initialize:
xi_nd_Po_s      = zeros(Nx1,Nx2);
xi_nd_Po_w      = zeros(Nx1,Nx2);
xi_nd_Po_e      = zeros(Nx1,Nx2);
xi_nd_Po_n      = zeros(Nx1,Nx2);

% Linear interpolate quantities at cell boundaries:
xi_nd_Po_s(i_in,j_in) = (xi_nd_Po_C(i_in,j_in) + xi_nd_Po_C(i_in,j_in-1))/2;  % south
xi_nd_Po_w(i_in,j_in) = (xi_nd_Po_C(i_in,j_in) + xi_nd_Po_C(i_in-1,j_in))/2;  % west
xi_nd_Po_e(i_in,j_in) = (xi_nd_Po_C(i_in+1,j_in) + xi_nd_Po_C(i_in,j_in))/2;  % east
xi_nd_Po_n(i_in,j_in) = (xi_nd_Po_C(i_in,j_in+1) + xi_nd_Po_C(i_in,j_in))/2;  % north 


% Couette term:
% Discretized with generic upwind scheme:
% Initialization:
xi_nd_Co_WW = zeros(Nx1,Nx2);
xi_nd_Co_W  = zeros(Nx1,Nx2);
xi_nd_Co_C  = zeros(Nx1,Nx2);
xi_nd_Co_E  = zeros(Nx1,Nx2);

if sld.flag_elastic 
    xi_nd_Co_h_WW = zeros(Nx1,Nx2);
    xi_nd_Co_h_W  = zeros(Nx1,Nx2);
    xi_nd_Co_h_C  = zeros(Nx1,Nx2);
    xi_nd_Co_h_E  = zeros(Nx1,Nx2);
end

% Definitions:
xi_nd_Co_cf                     = 12*ref.x1*ref.u*ref.mu/(ref.h^2*ref.p);
xi_nd_Co_W(i_nWW ,j_nWW)        = xi_nd_Co_cf .* (prop_nd.rho(i_nWW-1,j_nWW) .* prop_nd.u .* h_nd.h_m(i_nWW-1,j_nWW));
xi_nd_Co_C(i_nWW ,j_nWW)        = xi_nd_Co_cf .* (prop_nd.rho(i_nWW  ,j_nWW) .* prop_nd.u .* h_nd.h_m(i_nWW  ,j_nWW));
xi_nd_Co_E(i_nWW ,j_nWW)        = xi_nd_Co_cf .* (prop_nd.rho(i_nWW+1,j_nWW) .* prop_nd.u .* h_nd.h_m(i_nWW+1,j_nWW));

xi_nd_Co_WW(i_wWW,j_wWW)        = xi_nd_Co_cf .* (prop_nd.rho(i_wWW-2,j_wWW) .* prop_nd.u .* h_nd.h_m(i_wWW-2,j_wWW));
xi_nd_Co_W(i_wWW ,j_wWW)        = xi_nd_Co_cf .* (prop_nd.rho(i_wWW-1,j_wWW) .* prop_nd.u .* h_nd.h_m(i_wWW-1,j_wWW));
xi_nd_Co_C(i_wWW ,j_wWW)        = xi_nd_Co_cf .* (prop_nd.rho(i_wWW  ,j_wWW) .* prop_nd.u .* h_nd.h_m(i_wWW  ,j_wWW));
xi_nd_Co_E(i_wWW ,j_wWW)        = xi_nd_Co_cf .* (prop_nd.rho(i_wWW+1,j_wWW) .* prop_nd.u .* h_nd.h_m(i_wWW+1,j_wWW));

if sld.flag_elastic 
    xi_nd_Co_h_W(i_nWW ,j_nWW)  = xi_nd_Co_cf .* (prop_nd.rho(i_nWW-1,j_nWW) .* prop_nd.u);
    xi_nd_Co_h_C(i_nWW ,j_nWW)  = xi_nd_Co_cf .* (prop_nd.rho(i_nWW  ,j_nWW) .* prop_nd.u);
    xi_nd_Co_h_E(i_nWW ,j_nWW)  = xi_nd_Co_cf .* (prop_nd.rho(i_nWW+1,j_nWW) .* prop_nd.u);

    xi_nd_Co_h_WW(i_wWW,j_wWW)  = xi_nd_Co_cf .* (prop_nd.rho(i_wWW-2,j_wWW) .* prop_nd.u);
    xi_nd_Co_h_W(i_wWW ,j_wWW)  = xi_nd_Co_cf .* (prop_nd.rho(i_wWW-1,j_wWW) .* prop_nd.u);
    xi_nd_Co_h_C(i_wWW ,j_wWW)  = xi_nd_Co_cf .* (prop_nd.rho(i_wWW  ,j_wWW) .* prop_nd.u);
    xi_nd_Co_h_E(i_wWW ,j_wWW)  = xi_nd_Co_cf .* (prop_nd.rho(i_wWW+1,j_wWW) .* prop_nd.u);
end

% Generic upwind Couette term discretization is used for the Couette term (maximum order is 2 
% since midpoint rule is always used for integral flux approximations): 
% 1) 1.Order Upwind, 2) 2.Order QUICK, 3) 2.Order LUI, 4) 2.Order CUI
if alg.FBNS.flag_Couette_Discrete == 1
    a_Co = 0;
    b_Co = 1;
    c_Co = 0;
elseif alg.FBNS.flag_Couette_Discrete == 2
    a_Co = -1/8;
    b_Co =  6/8;
    c_Co =  3/8;
elseif alg.FBNS.flag_Couette_Discrete == 3
    a_Co = -1/2;
    b_Co =  3/2;
    c_Co =  0;
elseif alg.FBNS.flag_Couette_Discrete == 4
    a_Co = -1/6;
    b_Co =  5/6;
    c_Co =  2/6;
end

% Initialization:
xi_nd_Co_w  = zeros(Nx1,Nx2);
xi_nd_Co_e  = zeros(Nx1,Nx2);

% 1.Order upwind Couette discretization for west boundaries of cells inside domain but without a WestWest neighbour:
xi_nd_Co_w(i_nWW,j_nWW)         =      xi_nd_Co_W(   i_nWW,j_nWW);
% Generic upwind Couette discretization for other cell boundaries:
xi_nd_Co_e(i_nWW,j_nWW)         = a_Co*xi_nd_Co_W(   i_nWW,j_nWW);
xi_nd_Co_e(i_nWW,j_nWW)         = b_Co*xi_nd_Co_C(   i_nWW,j_nWW) + xi_nd_Co_e(  i_nWW,j_nWW);
xi_nd_Co_e(i_nWW,j_nWW)         = c_Co*xi_nd_Co_E(   i_nWW,j_nWW) + xi_nd_Co_e(  i_nWW,j_nWW);

xi_nd_Co_w(i_wWW,j_wWW)         = a_Co*xi_nd_Co_WW(  i_wWW,j_wWW);
xi_nd_Co_w(i_wWW,j_wWW)         = b_Co*xi_nd_Co_W(   i_wWW,j_wWW) + xi_nd_Co_w(  i_wWW,j_wWW);
xi_nd_Co_w(i_wWW,j_wWW)         = c_Co*xi_nd_Co_C(   i_wWW,j_wWW) + xi_nd_Co_w(  i_wWW,j_wWW);

xi_nd_Co_e(i_wWW,j_wWW)         = a_Co*xi_nd_Co_W(   i_wWW,j_wWW);
xi_nd_Co_e(i_wWW,j_wWW)         = b_Co*xi_nd_Co_C(   i_wWW,j_wWW) + xi_nd_Co_e(  i_wWW,j_wWW);
xi_nd_Co_e(i_wWW,j_wWW)         = c_Co*xi_nd_Co_E(   i_wWW,j_wWW) + xi_nd_Co_e(  i_wWW,j_wWW);



% Time dependent term:
% First timestep is the solution of the steady problem since 
% an initial solution in time is needed,
% while the other time steps consider the time dependent term through Euler
% implicit descretization.
if alg.flag_unsteady && (alg.time.it > 1)
       
    % Definitions:
    xi_nd_Ti_cf         = 12*ref.x1^2*ref.mu/(ref.t*ref.h^2*ref.p);
    xi_nd_Ti_C          = xi_nd_Ti_cf.*(     prop_nd.rho.*     h_nd.h_m);
    xi_nd_Ti_C_prev     = xi_nd_Ti_cf.*(prev.prop_nd.rho.*prev.h_nd.h_m);
    
    if sld.flag_elastic 
        xi_nd_Ti_h_C    = xi_nd_Ti_cf.*(prop_nd.rho);        
    end
end

% A:
% Initialize coefficients:
A_Po_cf_S      = zeros(Nx1,Nx2);
A_Po_cf_W      = zeros(Nx1,Nx2);
A_Po_cf_C      = zeros(Nx1,Nx2);
A_Po_cf_E      = zeros(Nx1,Nx2);
A_Po_cf_N      = zeros(Nx1,Nx2);

% Determine coefficients:        
A_Po_cf_S(i_in,j_in)  =  ((x_rat^2)* dx1_nd/dx2_nd)  .*xi_nd_Po_s(i_in,j_in);  % South 
A_Po_cf_W(i_in,j_in)  =  (           dx2_nd/dx1_nd)  .*xi_nd_Po_w(i_in,j_in);  % West
% Center :
A_Po_cf_C(i_in,j_in)  = -((x_rat^2)* dx1_nd/dx2_nd)  .*xi_nd_Po_s(i_in,j_in);
A_Po_cf_C(i_in,j_in)  = -(           dx2_nd/dx1_nd)  .*xi_nd_Po_w(i_in,j_in) + A_Po_cf_C(i_in,j_in);
A_Po_cf_C(i_in,j_in)  = -(           dx2_nd/dx1_nd)  .*xi_nd_Po_e(i_in,j_in) + A_Po_cf_C(i_in,j_in);
A_Po_cf_C(i_in,j_in)  = -((x_rat^2)* dx1_nd/dx2_nd)  .*xi_nd_Po_n(i_in,j_in) + A_Po_cf_C(i_in,j_in);

A_Po_cf_E(i_in,j_in)  =  (           dx2_nd/dx1_nd)  .*xi_nd_Po_e(i_in,j_in);  % East
A_Po_cf_N(i_in,j_in)  =  ((x_rat^2)* dx1_nd/dx2_nd)  .*xi_nd_Po_n(i_in,j_in);  % North 

% Initialize diagonals of A_Po:
% Notation: [line, column, value] of non-sparse matrix
A_Po_S     = zeros(N,3); % South diagonal
A_Po_W     = zeros(N,3); % West diagonal
A_Po_C     = zeros(N,3); % Center diagonal
A_Po_E     = zeros(N,3); % East diagonal
A_Po_N     = zeros(N,3); % North diagonal

% Fill diagonals of A with linear indices first along x1 then along x2:
% Notation: [line, column, value] of non-sparse matrix
% South:
A_Po_S(1:N_in,1)     = i_lin_in(:);            % South 
A_Po_S(1:N_in,2)     = i_lin_in(:) - Nx1;      % South 
A_Po_S(1:N_in,3)     = A_Po_cf_S(i_lin_in(:));    % South 
% West:
A_Po_W(1:N_in,1)     = i_lin_in(:);            % West 
A_Po_W(1:N_in,2)     = i_lin_in(:) - 1;        % West 
A_Po_W(1:N_in,3)     = A_Po_cf_W(i_lin_in(:));    % West
% Center :
A_Po_C(1:N_in,1)     = i_lin_in(:);            % Center 
A_Po_C(1:N_in,2)     = i_lin_in(:);            % Center
A_Po_C(1:N_in,3)     = A_Po_cf_C(i_lin_in(:));    % Center
% East:
A_Po_E(1:N_in,1)     = i_lin_in(:);            % East 
A_Po_E(1:N_in,2)     = i_lin_in(:) + 1;        % East 
A_Po_E(1:N_in,3)     = A_Po_cf_E(i_lin_in(:));    % East
% North:
A_Po_N(1:N_in,1)     = i_lin_in(:);            % North 
A_Po_N(1:N_in,2)     = i_lin_in(:) + Nx1;      % North
A_Po_N(1:N_in,3)     = A_Po_cf_N(i_lin_in(:));    % North

% B_Co:
% Initialize coefficients:
B_Co_cf_WW     = zeros(Nx1,Nx2);
B_Co_cf_W      = zeros(Nx1,Nx2);
B_Co_cf_C      = zeros(Nx1,Nx2);
B_Co_cf_E      = zeros(Nx1,Nx2);

% Determine coefficients:
% 1.Order upwind Couette discretization for west boundaries of cells inside domain but without a WestWest neighbour:
% West:
B_Co_cf_W( i_nWW,j_nWW)  = -(1   *dx2_nd).*xi_nd_Co_w(i_nWW,j_nWW);
B_Co_cf_W( i_nWW,j_nWW)  =  (a_Co*dx2_nd).*xi_nd_Co_e(i_nWW,j_nWW) + B_Co_cf_W( i_nWW,j_nWW);
% Center:
B_Co_cf_C( i_nWW,j_nWW)  =  (b_Co*dx2_nd).*xi_nd_Co_e(i_nWW,j_nWW);
% East
B_Co_cf_E( i_nWW,j_nWW)  =  (c_Co*dx2_nd).*xi_nd_Co_e(i_nWW,j_nWW);

% Generic upwind Couette discretization for other cell boundaries:
% WestWest:
B_Co_cf_WW(i_wWW,j_wWW)  = -(a_Co*dx2_nd).*xi_nd_Co_w(i_wWW,j_wWW);        
% West:
B_Co_cf_W( i_wWW,j_wWW)  = -(b_Co*dx2_nd).*xi_nd_Co_w(i_wWW,j_wWW);
B_Co_cf_W( i_wWW,j_wWW)  =  (a_Co*dx2_nd).*xi_nd_Co_e(i_wWW,j_wWW) + B_Co_cf_W( i_wWW,j_wWW);
% Center:
B_Co_cf_C( i_wWW,j_wWW)  = -(c_Co*dx2_nd).*xi_nd_Co_w(i_wWW,j_wWW);
B_Co_cf_C( i_wWW,j_wWW)  =  (b_Co*dx2_nd).*xi_nd_Co_e(i_wWW,j_wWW) + B_Co_cf_C( i_wWW,j_wWW);
% East
B_Co_cf_E( i_wWW,j_wWW)  =  (c_Co*dx2_nd).*xi_nd_Co_e(i_wWW,j_wWW);


% Initialize diagonals of B_Co:
B_Co_WW    = zeros(N_wWW,3); % WestWest diagonal
B_Co_W     = zeros(N_in ,3); % West diagonal
B_Co_C     = zeros(N_in ,3); % Center diagonal
B_Co_E     = zeros(N_in ,3); % East diagonal

% Define diagonals:
% WestWest:
B_Co_WW(1:N_wWW,1)     = i_lin_wWW(:);           
B_Co_WW(1:N_wWW,2)     = i_lin_wWW(:) - 2;      
B_Co_WW(1:N_wWW,3)     = B_Co_cf_WW(i_lin_wWW(:));   
% West:
B_Co_W( 1:N_in ,1)     = i_lin_in(:);         
B_Co_W( 1:N_in ,2)     = i_lin_in(:) - 1;     
B_Co_W( 1:N_in ,3)     = B_Co_cf_W( i_lin_in(:));  
% Center:
B_Co_C( 1:N_in ,1)     = i_lin_in(:);          
B_Co_C( 1:N_in ,2)     = i_lin_in(:);          
B_Co_C( 1:N_in ,3)     = B_Co_cf_C( i_lin_in(:));   
% Center:
B_Co_E( 1:N_in ,1)     = i_lin_in(:);          
B_Co_E( 1:N_in ,2)     = i_lin_in(:) + 1;          
B_Co_E( 1:N_in ,3)     = B_Co_cf_E( i_lin_in(:));  

% B_Ti:
if alg.flag_unsteady && (alg.time.it > 1)
    % Initialize:
    B_Ti_C     = zeros(N_in ,3); % Center diagonal
    % Define coefficient:
    B_Ti_cf_C = (dx1_nd*dx2_nd/alg.time.delta_t_nd).*xi_nd_Ti_C;
    % Fill:
    B_Ti_C( 1:N_in ,1)     = i_lin_in(:);          
    B_Ti_C( 1:N_in ,2)     = i_lin_in(:);          
    B_Ti_C( 1:N_in ,3)     = B_Ti_cf_C( i_lin_in(:));   
end

% c:
c = zeros(N,1);
% c_Co:
c(i_lin_in(:)) = dx2_nd*xi_nd_Co_w(i_lin_in(:)) - dx2_nd*xi_nd_Co_e(i_lin_in(:));
% c_Ti:
if alg.flag_unsteady && (alg.time.it > 1)
    c(i_lin_in(:)) = -(dx1_nd*dx2_nd/alg.time.delta_t_nd).*xi_nd_Ti_C(     i_lin_in(:))                                     + c(i_lin_in(:));
    c(i_lin_in(:)) =  (dx1_nd*dx2_nd/alg.time.delta_t_nd).*xi_nd_Ti_C_prev(i_lin_in(:)).*(1 - prev.sol.thet(i_lin_in(:)))   + c(i_lin_in(:));
end

% Boundaries with coefficients for better scaling of A:
% Off-main diagonal arrays are resized if they are finished
% South boundary:
k_l                             = N_in + 1;
k_u                             = N_in + Nx1;
if fld_nd.flag_p_bound == 0
    % Dirichlet:
    A_Po_C(k_l:k_u,1)              = i_lin(1:Nx1,1);
    A_Po_C(k_l:k_u,2)              = i_lin(1:Nx1,1);
    A_Po_C(k_l:k_u,3)              =  ((x_rat^2)* dx1_nd/dx2_nd).* xi_nd_Po_C(1:Nx1,1);  
    
    c(i_lin(1:Nx1,1))           =  ((x_rat^2)* dx1_nd/dx2_nd).* xi_nd_Po_C(1:Nx1,1);
    c(i_lin(1:Nx1,1))           = -(fld_nd.p_S_side - fld_nd.p_cav) .*c(i_lin(1:Nx1,1));
    
    A_Po_N = A_Po_N(1:N_in,1:3);       
    
elseif fld_nd.flag_p_bound == 1 || fld_nd.flag_p_bound == 2
    % Neumann:
    A_Po_C(k_l:k_u,1)              = i_lin(1:Nx1,1);
    A_Po_C(k_l:k_u,2)              = i_lin(1:Nx1,1);
    A_Po_C(k_l:k_u,3)              =  ((x_rat^2)* dx1_nd/dx2_nd).* xi_nd_Po_C(1:Nx1,1);
    
    A_Po_N(N_in+1:N_in+Nx1,1)      = i_lin(1:Nx1,1);
    A_Po_N(N_in+1:N_in+Nx1,2)      = i_lin(1:Nx1,1) + Nx1;
    A_Po_N(N_in+1:N_in+Nx1,3)      = -((x_rat^2)* dx1_nd/dx2_nd).* xi_nd_Po_C(1:Nx1,1);
    A_Po_N                         = A_Po_N(1:N_in + Nx1,1:3);
    
end

% West boundary:
k_l                             = k_u + 1;
k_u                             = k_u + Nx2 - 2;
% Dirichlet:
A_Po_C(k_l:k_u,1)                  = i_lin(1,2:Nx2-1)';
A_Po_C(k_l:k_u,2)                  = i_lin(1,2:Nx2-1)';    
A_Po_C(k_l:k_u,3)                  =  (            dx2_nd/dx1_nd).*(xi_nd_Po_C(1,2:Nx2-1)');  

c(i_lin(1,2:Nx2-1)')            =  (            dx2_nd/dx1_nd).*(xi_nd_Po_C(1,2:Nx2-1)');
c(i_lin(1,2:Nx2-1)')            = -(fld_nd.p_W_side - fld_nd.p_cav) .*c(i_lin(1,2:Nx2-1)');

A_Po_E                             = A_Po_E(1:N_in,1:3);

% East boundary:
k_l                             = k_u + 1;
k_u                             = k_u + Nx2 - 2;
if fld_nd.flag_p_bound == 0 || fld_nd.flag_p_bound == 1
    % Dirichlet:
    A_Po_C(k_l:k_u,1)              = i_lin(Nx1,2:Nx2-1)';
    A_Po_C(k_l:k_u,2)              = i_lin(Nx1,2:Nx2-1)';
    A_Po_C(k_l:k_u,3)              =  (            dx2_nd/dx1_nd).*(xi_nd_Po_C(Nx1,2:Nx2-1)');  
    
    c(i_lin(Nx1,2:Nx2-1)')      =  (            dx2_nd/dx1_nd).*(xi_nd_Po_C(Nx1,2:Nx2-1)');
    c(i_lin(Nx1,2:Nx2-1)')      = -(fld_nd.p_E_side - fld_nd.p_cav) .*c(i_lin(Nx1,2:Nx2-1)');
    
    A_Po_W                         = A_Po_W(1:N_in,1:3);
    
elseif fld_nd.flag_p_bound == 2
    % Neumann:
    A_Po_C(k_l:k_u,1)              = i_lin(Nx1,2:Nx2-1)';
    A_Po_C(k_l:k_u,2)              = i_lin(Nx1,2:Nx2-1)';
    A_Po_C(k_l:k_u,3)              =  (            dx2_nd/dx1_nd).*(xi_nd_Po_C(Nx1,2:Nx2-1)');  
    
    A_Po_W(N_in+1:N_in+Nx2-2,1)    = i_lin(Nx1,2:Nx2-1)';
    A_Po_W(N_in+1:N_in+Nx2-2,2)    = i_lin(Nx1,2:Nx2-1)' - 1;
    A_Po_W(N_in+1:N_in+Nx2-2,3)    = -(            dx2_nd/dx1_nd).*(xi_nd_Po_C(Nx1,2:Nx2-1)');
    A_Po_W                         = A_Po_W(1:N_in+Nx2-2,1:3);
    
end

% North boundary:
k_l                             = k_u + 1;
k_u                             = k_u + Nx1;
if fld_nd.flag_p_bound == 0
    % Dirichlet:
    A_Po_C(k_l:k_u,1)              = i_lin(1:Nx1,Nx2);
    A_Po_C(k_l:k_u,2)              = i_lin(1:Nx1,Nx2);
    A_Po_C(k_l:k_u,3)              =  ((x_rat^2)* dx1_nd/dx2_nd).* xi_nd_Po_C(1:Nx1,Nx2); 
    
    c(i_lin(1:Nx1,Nx2))         =  ((x_rat^2)* dx1_nd/dx2_nd).* xi_nd_Po_C(1:Nx1,Nx2);
    c(i_lin(1:Nx1,Nx2))         = -(fld_nd.p_N_side - fld_nd.p_cav) .*c(i_lin(1:Nx1,Nx2));
    
    A_Po_S                         = A_Po_S(1:N_in,1:3);  
    
elseif fld_nd.flag_p_bound == 1 || fld_nd.flag_p_bound == 2
    % Neumann:
    A_Po_C(k_l:k_u,1)              = i_lin(1:Nx1,Nx2);
    A_Po_C(k_l:k_u,2)              = i_lin(1:Nx1,Nx2);
    A_Po_C(k_l:k_u,3)              =  ((x_rat^2)* dx1_nd/dx2_nd).* xi_nd_Po_C(1:Nx1,Nx2);
    
    A_Po_S(N_in+1:N_in+Nx1,1)      = i_lin(1:Nx1,Nx2);
    A_Po_S(N_in+1:N_in+Nx1,2)      = i_lin(1:Nx1,Nx2) - Nx1;
    A_Po_S(N_in+1:N_in+Nx1,3)      = -((x_rat^2)* dx1_nd/dx2_nd).* xi_nd_Po_C(1:Nx1,Nx2);
    A_Po_S                         = A_Po_S(1:N_in + Nx1,1:3);
    
end

% Assemble A in sparse notation:
A = sparse(A_Po_S(:,1),A_Po_S(:,2),A_Po_S(:,3),N,N);
A = sparse(A_Po_W(:,1),A_Po_W(:,2),A_Po_W(:,3),N,N) + A; 
A = sparse(A_Po_C(:,1),A_Po_C(:,2),A_Po_C(:,3),N,N) + A;
A = sparse(A_Po_E(:,1),A_Po_E(:,2),A_Po_E(:,3),N,N) + A;
A = sparse(A_Po_N(:,1),A_Po_N(:,2),A_Po_N(:,3),N,N) + A;

% Assemble B in sparse notation:
B = sparse(B_Co_WW(:,1),B_Co_WW(:,2),B_Co_WW(:,3),N,N);
B = sparse(B_Co_W( :,1),B_Co_W( :,2),B_Co_W( :,3),N,N) + B;
B = sparse(B_Co_C( :,1),B_Co_C( :,2),B_Co_C( :,3),N,N) + B;
B = sparse(B_Co_E( :,1),B_Co_E( :,2),B_Co_E( :,3),N,N) + B;
if alg.flag_unsteady && (alg.time.it > 1)
    B = sparse(B_Ti_C( :,1),B_Ti_C( :,2),B_Ti_C( :,3),N,N) + B;
end


% Construct J_Gp_nd:
J_Gp_nd = A;

if sld.flag_elastic         
    % Construction of A_h with Kernel entries: 
    % Construct sparse A_h containing the five main diagonals:
    % West:
    A_Co_h_W       = construct_A_Co_h_entries(-1, 0); 
    A_Co_h         = sparse(A_Co_h_W(:,1),A_Co_h_W(:,2),A_Co_h_W(:,3),N,N);
    % Center:
    A_Co_h_C       = construct_A_Co_h_entries( 0, 0);
    A_Co_h         = sparse(A_Co_h_C(:,1),A_Co_h_C(:,2),A_Co_h_C(:,3),N,N) + A_Co_h;
    % South:
    if alg.FBNS.flag_A_h_S
        A_Co_h_S   = construct_A_Co_h_entries( 0,-1); 
        A_Co_h     = sparse(A_Co_h_S(:,1),A_Co_h_S(:,2),A_Co_h_S(:,3),N,N)  + A_Co_h;
    end
    % WestWest:
    if alg.FBNS.flag_A_h_WW
        A_Co_h_WW  = construct_A_Co_h_entries(-2, 0);
        A_Co_h     = sparse(A_Co_h_WW(:,1),A_Co_h_WW(:,2),A_Co_h_WW(:,3),N,N) + A_Co_h;
    end
    % East:   
    if alg.FBNS.flag_A_h_E
        A_Co_h_E   = construct_A_Co_h_entries( 1, 0); 
        A_Co_h     = sparse(A_Co_h_E(:,1),A_Co_h_E(:,2),A_Co_h_E(:,3),N,N) + A_Co_h;
    end
    % North:
    if alg.FBNS.flag_A_h_N
        A_Co_h_N   = construct_A_Co_h_entries( 0, 1); 
        A_Co_h     = sparse(A_Co_h_N(:,1),A_Co_h_N(:,2),A_Co_h_N(:,3),N,N) + A_Co_h;
    end
    
    % Construct J_Gp:
    J_Gp_nd = J_Gp_nd + A_Co_h;
    
    
    % Unsteady term:
    if alg.flag_unsteady && (alg.time.it > 1)        
        % The kernel is saved in a shifted way such that the kernel Center is at entry (1,1).
        % This means that e.g. the entry to the West of the Center is saved at (2*Nx1,1).
        % To extract the required kernel entries, the modulo operator can be used 
        % to identify the corresponding indices:
        i_W  = mod(-1, 2*Nx1) + 1;
        i_C  = mod( 0, 2*Nx1) + 1;
        i_E  = mod( 1, 2*Nx1) + 1;
        i_EE = mod( 2, 2*Nx1) + 1;

        j_S  = mod(-1, 2*Nx2) + 1;
        j_C  = mod( 0, 2*Nx2) + 1;
        j_N  = mod( 1, 2*Nx2) + 1;

        % Extract Kernel entries:
        K_nd_S     = h_nd.Kernel(i_C ,j_S);
        K_nd_W     = h_nd.Kernel(i_W ,j_C);
        K_nd_C     = h_nd.Kernel(i_C ,j_C);
        K_nd_E     = h_nd.Kernel(i_E ,j_C);
        K_nd_EE    = h_nd.Kernel(i_EE,j_C);
        K_nd_N     = h_nd.Kernel(i_C ,j_N);

        % Fill diagonals of A_h_Ti with linear indices first along x1 then along x2:
        % Notation: [line, column, value] of non-sparse matrix
        % West:
        A_Ti_h_W = zeros(N_in,3);
        A_Ti_h_W(1:N_in,1) = i_lin_in(:);
        A_Ti_h_W(1:N_in,2) = i_lin_in(:) - 1;
        A_Ti_h_W(1:N_in,3) = -(dx1_nd*dx2_nd/alg.time.delta_t_nd*K_nd_E ).*xi_nd_Ti_h_C( i_lin_in(:));
        A_Ti_h          = sparse(A_Ti_h_W(:,1),A_Ti_h_W(:,2),A_Ti_h_W(:,3),N,N);
        % Center:
        A_Ti_h_C = zeros(N_in,3);
        A_Ti_h_C(1:N_in,1) = i_lin_in(:);
        A_Ti_h_C(1:N_in,2) = i_lin_in(:);
        A_Ti_h_C(1:N_in,3) = -(dx1_nd*dx2_nd/alg.time.delta_t_nd*K_nd_C ).*xi_nd_Ti_h_C( i_lin_in(:));
        A_Ti_h          = sparse(A_Ti_h_C(:,1),A_Ti_h_C(:,2),A_Ti_h_C(:,3),N,N) + A_Ti_h;
        % South:
        if alg.FBNS.flag_A_h_S
            A_Ti_h_S = zeros(N_in,3);
            A_Ti_h_S(1:N_in,1)  = i_lin_in(:);
            A_Ti_h_S(1:N_in,2)  = i_lin_in(:) - Nx1;
            A_Ti_h_S(1:N_in,3)  = -(dx1_nd*dx2_nd/alg.time.delta_t_nd*K_nd_N ).*xi_nd_Ti_h_C( i_lin_in(:));
            A_Ti_h              = sparse(A_Ti_h_S( :,1),A_Ti_h_S( :,2),A_Ti_h_S( :,3),N,N) + A_Ti_h;
        end
        % WestWest:
        if alg.FBNS.flag_A_h_WW
            A_Ti_h_WW = zeros(N_in,3);
            A_Ti_h_WW(1:N_in,1) = i_lin_in(:);
            A_Ti_h_WW(1:N_in,2) = i_lin_in(:) - 2;
            A_Ti_h_WW(1:N_in,3) = -(dx1_nd*dx2_nd/alg.time.delta_t_nd*K_nd_EE).*xi_nd_Ti_h_C( i_lin_in(:));
            A_Ti_h              = sparse(A_Ti_h_WW(:,1),A_Ti_h_WW(:,2),A_Ti_h_WW(:,3),N,N) + A_Ti_h;
       end
        
        if alg.FBNS.flag_A_h_E
            A_Ti_h_E = zeros(N_in,3);
            A_Ti_h_E(1:N_in,1)  = i_lin_in(:);
            A_Ti_h_E(1:N_in,2)  = i_lin_in(:) + 1;
            A_Ti_h_E(1:N_in,3)  = -(dx1_nd*dx2_nd/alg.time.delta_t_nd*K_nd_W ).*xi_nd_Ti_h_C( i_lin_in(:));
            A_Ti_h              = sparse(A_Ti_h_E( :,1),A_Ti_h_E( :,2),A_Ti_h_E( :,3),N,N) + A_Ti_h;
       end
        if alg.FBNS.flag_A_h_N
            A_Ti_h_N = zeros(N_in,3);
            A_Ti_h_N(1:N_in,1)  = i_lin_in(:);
            A_Ti_h_N(1:N_in,2)  = i_lin_in(:) + Nx1;
            A_Ti_h_N(1:N_in,3)  = -(dx1_nd*dx2_nd/alg.time.delta_t_nd*K_nd_S ).*xi_nd_Ti_h_C( i_lin_in(:));
            A_Ti_h              = sparse(A_Ti_h_N( :,1),A_Ti_h_N( :,2),A_Ti_h_N( :,3),N,N) + A_Ti_h;
        end
        
    % Construct J_Gp:
    J_Gp_nd = J_Gp_nd + A_Ti_h;
    end
end

    function A_Co_h_gen = construct_A_Co_h_entries(i_d,j_d)   
    % The kernel is saved in a shifted way such that the kernel Center is at entry (1,1).
    % This means that e.g. the entry to the West of the Center is saved at (2*Nx1,1).
    % To extract the required kernel entries, the modulo operator can be used 
    % to identify the corresponding indices:
    i_WW = mod(-i_d - 2, 2*Nx1) + 1;
    i_W  = mod(-i_d - 1, 2*Nx1) + 1;
    i_C  = mod(-i_d    , 2*Nx1) + 1;
    i_E  = mod(-i_d + 1, 2*Nx1) + 1;
    
    j_C  = mod(-j_d    , 2*Nx2) + 1;

    % Extract Kernel entries:
    K_nd_WW    = h_nd.Kernel(i_WW,j_C);
    K_nd_W     = h_nd.Kernel(i_W ,j_C);
    K_nd_C     = h_nd.Kernel(i_C ,j_C);
    K_nd_E     = h_nd.Kernel(i_E ,j_C);
    
        
    % Fill diagonals of A_h with linear indices first along x1 then along x2:
    % Notation: [line, column, value] of non-sparse matrix
    A_Co_h_gen = zeros(N_in,3);
    
    % 1.Order upwind Couette discretization for west boundaries of cells inside domain but without a WestWest neighbour:
    % South:
    A_Co_h_gen(1:N_nWW,1) = i_lin_nWW(:);
    A_Co_h_gen(1:N_nWW,2) = i_lin_nWW(:) + i_d + Nx1*j_d;
    A_Co_h_gen(1:N_nWW,3) =  (1   *dx2_nd*K_nd_W ).*xi_nd_Co_h_W( i_lin_nWW(:));
    A_Co_h_gen(1:N_nWW,3) = -(a_Co*dx2_nd*K_nd_W ).*xi_nd_Co_h_W( i_lin_nWW(:)) + A_Co_h_gen(1:N_nWW,3);
    A_Co_h_gen(1:N_nWW,3) = -(b_Co*dx2_nd*K_nd_C ).*xi_nd_Co_h_C( i_lin_nWW(:)) + A_Co_h_gen(1:N_nWW,3);
    A_Co_h_gen(1:N_nWW,3) = -(c_Co*dx2_nd*K_nd_E ).*xi_nd_Co_h_E( i_lin_nWW(:)) + A_Co_h_gen(1:N_nWW,3);
    
    % Generic upwind Couette discretization for other cell boundaries:
    li = 1 + N_nWW;
    % South:
    A_Co_h_gen(li:N_in,1) = i_lin_wWW(:);
    A_Co_h_gen(li:N_in,2) = i_lin_wWW(:) + i_d + Nx1*j_d;
    A_Co_h_gen(li:N_in,3) =  (a_Co*dx2_nd*K_nd_WW).*xi_nd_Co_h_WW(i_lin_wWW(:));
    A_Co_h_gen(li:N_in,3) =  (b_Co*dx2_nd*K_nd_W ).*xi_nd_Co_h_W( i_lin_wWW(:)) + A_Co_h_gen(li:N_in,3);
    A_Co_h_gen(li:N_in,3) =  (c_Co*dx2_nd*K_nd_C ).*xi_nd_Co_h_C( i_lin_wWW(:)) + A_Co_h_gen(li:N_in,3);
    A_Co_h_gen(li:N_in,3) = -(a_Co*dx2_nd*K_nd_W ).*xi_nd_Co_h_W( i_lin_wWW(:)) + A_Co_h_gen(li:N_in,3);
    A_Co_h_gen(li:N_in,3) = -(b_Co*dx2_nd*K_nd_C ).*xi_nd_Co_h_C( i_lin_wWW(:)) + A_Co_h_gen(li:N_in,3);
    A_Co_h_gen(li:N_in,3) = -(c_Co*dx2_nd*K_nd_E ).*xi_nd_Co_h_E( i_lin_wWW(:)) + A_Co_h_gen(li:N_in,3);
        
    end
end



function [G,F,J_Fp_C,J_Ft_C] = FBNS_compute_res_nd(A,B,c,p_red_nd,thet_red,h_nd,fld_nd)
% This function computes the residual of the Reynolds equation G = A*p_nd_red + B*thet_red + c,
% the residual of the Fischer-Burmeister equation F = p_nd_red + thet_red - sqrt(p_nd_red^2 + thet_red^2)
% and the center lines of the Jacobians J_Fp_nd_C and J_Fthet_C.
% The boundary conditions of thet_red are incorporated in F, J_Fp_nd_C and J_Fthet_C and the boundary
% conditions of p_nd_red are already incorporated in G, A, B and c
 
% Shorter definitions:
N       = h_nd.N;
Nx1     = h_nd.Nx1;
Nx2     = h_nd.Nx2;
i_lin   = h_nd.i_lin;

% Inside domain excluding boundaries:
i           = 2:Nx1-1;           % regular indices in x1 direction
j           = 2:Nx2-1;           % regular indices in x2 direction
i_lin_in    = i_lin(i,j);        % linear indices

% Residual of Reynolds equation:
G = A*p_red_nd + B*thet_red + c;

% Residual of cavitation condition:
F = zeros(N,1);
F(i_lin_in(:)) =                         p_red_nd(i_lin_in(:))     +  thet_red(i_lin_in(:));
F(i_lin_in(:)) =  F(i_lin_in(:)) - sqrt((p_red_nd(i_lin_in(:))).^2 + (thet_red(i_lin_in(:))).^2);

% Boundaries:
% Boundary condition for p_red_nd is already considered in G through A, B and c
% Boundary conditions are put up in a way that center diagnoal entries 
% of J_Fp_nd and J_Fthet are >=0 or else the swapping process will not work 
% South boundary:
if fld_nd.flag_thet_bound == 0
    % Dirichlet:
    F(i_lin(1:Nx1,1))       =   thet_red(i_lin(1:Nx1,1)); 
    F(i_lin(1:Nx1,1))       = - fld_nd.thet_S_side                  + F(i_lin(1:Nx1,1)); 
    
elseif fld_nd.flag_thet_bound == 1 || fld_nd.flag_thet_bound == 2
    % Neumann:
    F(i_lin(1:Nx1,1))       =   thet_red(i_lin(1:Nx1,1));
    F(i_lin(1:Nx1,1))       = - thet_red(i_lin(1:Nx1,2))            + F(i_lin(1:Nx1,1));
    
end

% West boundary:
% Dirichlet:
F(i_lin(1,2:Nx2-1)')        =   thet_red(i_lin(1,2:Nx2-1)');
F(i_lin(1,2:Nx2-1)')        = - fld_nd.thet_W_side                  + F(i_lin(1,2:Nx2-1)');

% East boundary:
if fld_nd.flag_thet_bound == 0 || fld_nd.flag_thet_bound == 1
    % Dirichlet:
    F(i_lin(Nx1,2:Nx2-1)')  =   thet_red(i_lin(Nx1,2:Nx2-1)');
    F(i_lin(Nx1,2:Nx2-1)')  = - fld_nd.thet_E_side                  + F(i_lin(Nx1,2:Nx2-1)');
    
elseif fld_nd.flag_thet_bound == 2
    % Neumann:
    F(i_lin(Nx1,2:Nx2-1)')  =   thet_red(i_lin(Nx1,2:Nx2-1)');
    F(i_lin(Nx1,2:Nx2-1)')  = - thet_red(i_lin(Nx1-1,2:Nx2-1)')     + F(i_lin(Nx1,2:Nx2-1)');
    
end

% North boundary:
if fld_nd.flag_thet_bound == 0
    % Dirichlet:
    F(i_lin(1:Nx1,Nx2))     =   thet_red(i_lin(1:Nx1,Nx2));
    F(i_lin(1:Nx1,Nx2))     = - fld_nd.thet_N_side                  + F(i_lin(1:Nx1,Nx2));
    
elseif fld_nd.flag_thet_bound == 1 || fld_nd.flag_thet_bound == 2
    % Neumann:
    F(i_lin(1:Nx1,Nx2))     =   thet_red(i_lin(1:Nx1,Nx2));
    F(i_lin(1:Nx1,Nx2))     = - thet_red(i_lin(1:Nx1,Nx2-1))        + F(i_lin(1:Nx1,Nx2));
    
end

% J_Fp_nd and J_Fthet are approximated by their center diagonals of length N
% At the Theta Neumann boundaries, J_Fthet is reduced to the center diagonal
% All entries of J_Fp_nd and J_Fthet must be >=0 or else the swapping process will not work 
% Jacobian approximation of cavitation condition:
J_Fp_C = zeros(N,1);
J_Ft_C = zeros(N,1);
% FBNS algorithm fails if p_red_nd = thet_red = 0 in J_F, to rectify this:
p_red_nd_aux = p_red_nd;
p_red_nd_aux((p_red_nd< eps) & (p_red_nd>=0)) =  eps;
p_red_nd_aux((p_red_nd>-eps) & (p_red_nd< 0)) = -eps;
thet_red_aux = thet_red;
thet_red_aux((thet_red< eps) & (thet_red>=0)) =  eps;
thet_red_aux((thet_red>-eps) & (thet_red< 0)) = -eps;
J_Fp_C(i_lin_in(:)) = 1 - (p_red_nd_aux(i_lin_in(:)))                                       ./ ...
                     sqrt((p_red_nd_aux(i_lin_in(:))).^2 + (thet_red_aux(i_lin_in(:))).^2);
J_Ft_C(i_lin_in(:)) = 1 -                                  (thet_red_aux(i_lin_in(:)))      ./ ...
                     sqrt((p_red_nd_aux(i_lin_in(:))).^2 + (thet_red_aux(i_lin_in(:))).^2);

% Boundaries:
% Boundary condition for p_red_nd is already considered in G through A, B and c
% Remember that J_F is approximated only by its main diagonal, even though
% Neumann Boundaries would actually result in more entries
% South boundary:
J_Ft_C(i_lin(1:Nx1,1))          = 1;
% West boundary:
J_Ft_C(i_lin(1,2:Nx2-1)')       = 1;
% East boundary:
J_Ft_C(i_lin(Nx1,2:Nx2-1)')     = 1;
% North boundary:
J_Ft_C(i_lin(1:Nx1,Nx2))        = 1;
end

function [A_G,B_G,A_F,B_F,sw_cond] = FBNS_swapping(J_Gp_nd,J_Gthet,J_Fp_nd_C,J_Fthet_C,h_nd)
% Shorter definitions:
Nx1 = h_nd.Nx1;
N   = h_nd.N;
i   = 1:h_nd.N;                    % [-]   counter for rows
j   = 1:h_nd.N;                    % [-]   counter for columns

% Extract diagonals:    
% Extract diagonals of J_Gp_nd and J_Gthet:
% spdgiags extracts the diagonal of a given matrix. The index of the
% resulting vector corresponds to the column index of the  input matrix
J_Gp_nd_S  = spdiags(J_Gp_nd,-Nx1);         % South diagonal
J_Gp_nd_WW = spdiags(J_Gp_nd,-2);           % West diagonal
J_Gp_nd_W  = spdiags(J_Gp_nd,-1);           % West diagonal
J_Gp_nd_C  = spdiags(J_Gp_nd,0);            % Center diagonal
J_Gp_nd_E  = spdiags(J_Gp_nd,1);            % East diagonal
J_Gp_nd_N  = spdiags(J_Gp_nd,Nx1);          % North diagonal

J_Gthet_WW = spdiags(J_Gthet,-2);           % WestWest diagonal
J_Gthet_W  = spdiags(J_Gthet,-1);           % West diagonal
J_Gthet_C  = spdiags(J_Gthet,0);            % Center diagonal
J_Gthet_E  = spdiags(J_Gthet,1);            % Center diagonal

% Construct A_F, B_F, A_G and B_G:
% Initialize diagonals:
A_G_S   = J_Gp_nd_S;                        % South diagonal
A_G_WW  = J_Gp_nd_WW;                       % WestWest diagonal
A_G_W   = J_Gp_nd_W;                        % West diagonal
A_G_C   = J_Gp_nd_C;                        % Center diagonal
A_G_E   = J_Gp_nd_E;                        % East diagonal
A_G_N   = J_Gp_nd_N;                        % North diagonal

B_G_S   = zeros(N,1);                       % South diagonal, 0 since J_Gthet does not have a South diagonal
B_G_WW  = J_Gthet_WW;                       % West diagonal
B_G_W   = J_Gthet_W;                        % West diagonal
B_G_C   = J_Gthet_C;                        % Center diagonal
B_G_E   = J_Gthet_E;                        % East diagonal
B_G_N   = zeros(N,1);                       % North diagonal, 0 since J_Gthet does not have a North diagonal

A_F_C   = J_Fp_nd_C;                        % Center diagonal
B_F_C   = J_Fthet_C;                        % Center diagonal

% Swap condition to ensure that A_F is invertible and well conditioned:
sw_cond = J_Fp_nd_C<J_Fthet_C;       

% Swap columns of matrices:
% Swap center diagonals:
A_G_C(sw_cond)  = J_Gthet_C(sw_cond);
B_G_C(sw_cond)  = J_Gp_nd_C(sw_cond);
A_F_C(sw_cond)  = J_Fthet_C(sw_cond);
B_F_C(sw_cond)  = J_Fp_nd_C(sw_cond);

% Swap other diagonals:
A_G_S(sw_cond)  = 0;       % because the B_G_S zero before swapping
B_G_S(sw_cond)  = J_Gp_nd_S(sw_cond);
A_G_WW(sw_cond) = J_Gthet_WW(sw_cond);
B_G_WW(sw_cond) = J_Gp_nd_WW(sw_cond);
A_G_W(sw_cond)  = J_Gthet_W(sw_cond);
B_G_W(sw_cond)  = J_Gp_nd_W(sw_cond);
A_G_E(sw_cond)  = J_Gthet_E(sw_cond);
B_G_E(sw_cond)  = J_Gp_nd_E(sw_cond);
A_G_N(sw_cond)  = 0;
B_G_N(sw_cond)  = J_Gp_nd_N(sw_cond);

% Construct sparse matrizes:
A_G =   sparse(i(1+Nx1:N)   ,j(1:N-Nx1) ,A_G_S(1:N-Nx1) ,N,N);
A_G =   sparse(i(3:N)       ,j(1:N-2)   ,A_G_WW(1:N-2)  ,N,N) + A_G;
A_G =   sparse(i(2:N)       ,j(1:N-1)   ,A_G_W(1:N-1)   ,N,N) + A_G;
A_G =   sparse(i            ,j          ,A_G_C          ,N,N) + A_G;
A_G =   sparse(i(1:N-1)     ,j(2:N)     ,A_G_E(2:N)     ,N,N) + A_G;
A_G =   sparse(i(1:N-Nx1)   ,j(1+Nx1:N) ,A_G_N(1+Nx1:N) ,N,N) + A_G;
      
B_G =   sparse(i(1+Nx1:N)   ,j(1:N-Nx1) ,B_G_S(1:N-Nx1) ,N,N); 
B_G =   sparse(i(3:N)       ,j(1:N-2)   ,B_G_WW(1:N-2)  ,N,N) + B_G;
B_G =   sparse(i(2:N)       ,j(1:N-1)   ,B_G_W(1:N-1)   ,N,N) + B_G;
B_G =   sparse(i            ,j          ,B_G_C          ,N,N) + B_G;
B_G =   sparse(i(1:N-1)     ,j(2:N)     ,B_G_E(2:N)     ,N,N) + B_G;
B_G =   sparse(i(1:N-Nx1)   ,j(1+Nx1:N) ,B_G_N(1+Nx1:N) ,N,N) + B_G;

A_F =   sparse(i            ,j          ,A_F_C          ,N,N);
B_F =   sparse(i            ,j          ,B_F_C          ,N,N);
end

function [G_old,F_old,res] = FBNS_determine_residual(res,alg,delta_p_nd,delta_thet,G,F,G_old,F_old,sol,ref,opc)
% Determine different residuals:
res.FBNS.delta_p_nd_mean(alg.it_tot,1)   = mean(abs(delta_p_nd));
res.FBNS.delta_thet_mean(alg.it_tot,1)   = mean(abs(delta_thet));
res.FBNS.delta_p_nd_max(alg.it_tot,1)    = max(abs(delta_p_nd));
res.FBNS.delta_thet_max(alg.it_tot,1)    = max(abs(delta_thet));
res.FBNS.delta_G_max(alg.it_tot,1)       = max(abs(G - G_old));
res.FBNS.delta_F_max(alg.it_tot,1)       = max(abs(F - F_old));
res.FBNS.G_max(alg.it_tot,1)             = max(abs(G));
res.FBNS.F_max(alg.it_tot,1)             = max(abs(F));

% Save for next iteration:
G_old = G;
F_old = F;
    
% Determine which resiudual is used to evaluate FBNS algorithm:
if alg.FBNS.flag_res == 0
    res_aux     = zeros(2,1);
    res_aux(1)  = res.FBNS.delta_p_nd_mean(alg.it_tot,1);
    res_aux(2)  = res.FBNS.delta_thet_mean(alg.it_tot,1);
    % Find maximum:
    res.FBNS.FBNS(alg.it_tot,1) = max(res_aux(:));
    
elseif alg.FBNS.flag_res == 1
    res_aux     = zeros(2,1);
    res_aux(1)  = res.FBNS.delta_p_nd_max(alg.it_tot,1);
    res_aux(2)  = res.FBNS.delta_thet_max(alg.it_tot,1);
    % Find maximum:
    res.FBNS.FBNS(alg.it_tot,1) = max(res_aux(:));
    
elseif alg.FBNS.flag_res == 3
    res.FBNS.FBNS(alg.it_tot,1) = res.FBNS.delta_G_max(alg.it_tot,1);
    
elseif alg.FBNS.flag_res == 4
    res.FBNS.FBNS(alg.it_tot,1) = res.FBNS.G_max(alg.it_tot,1);
    
elseif alg.FBNS.flag_res == 5
    res_aux     = zeros(2,1);
    res_aux(1)  = res.FBNS.delta_G_max(alg.it_tot,1);
    res_aux(2)  = res.FBNS.G_max(alg.it_tot,1);
    % Find maximum:
    res.FBNS.FBNS(alg.it_tot,1) = max(res_aux(:));
    
elseif alg.FBNS.flag_res == 6
    res.FBNS.FBNS(alg.it_tot,1) = res.FBNS.delta_p_nd_mean(alg.it_tot,1);
    
elseif alg.FBNS.flag_res == 7
    res.FBNS.FBNS(alg.it_tot,1) = res.FBNS.delta_p_nd_max(alg.it_tot,1);
    
elseif alg.FBNS.flag_res == 8
    res_aux     = zeros(3,1);
    res_aux(1)  = res.FBNS.G_max(alg.it_tot,1);
    res_aux(2)  = res.FBNS.delta_p_nd_mean(alg.it_tot,1);
    res_aux(3)  = res.FBNS.delta_thet_mean(alg.it_tot,1);
    % Find maximum:
    res.FBNS.FBNS(alg.it_tot,1) = max(res_aux(:));
    
elseif alg.FBNS.flag_res == 9
    res_aux     = zeros(3,1);
    res_aux(1)  = res.FBNS.G_max(alg.it_tot,1) ;
    res_aux(2)  = res.FBNS.F_max(alg.it_tot,1);
    % Find maximum:
    res.FBNS.FBNS(alg.it_tot,1) = max(res_aux(:));
    
elseif alg.FBNS.flag_res == 10
    res_aux     = zeros(6,1);
    res_aux(1)  = res.FBNS.G_max(alg.it_tot,1);
    res_aux(2)  = res.FBNS.F_max(alg.it_tot,1);
    res_aux(3)  = res.FBNS.delta_G_max(alg.it_tot,1);
    res_aux(4)  = res.FBNS.delta_F_max(alg.it_tot,1);
    res_aux(5)  = res.FBNS.delta_p_nd_max(alg.it_tot,1);
    res_aux(6)  = res.FBNS.delta_thet_max(alg.it_tot,1);
    % Find maximum:
    res.FBNS.FBNS(alg.it_tot,1) = max(res_aux(:));
    
elseif alg.FBNS.flag_res == 11
    res_aux     = zeros(3,1);
    res_aux(1)  = res.FBNS.G_max(alg.it_tot,1);
    res_aux(2)  = res.FBNS.delta_G_max(alg.it_tot,1);
    res_aux(3)  = res.FBNS.delta_p_nd_max(alg.it_tot,1);
    % Find maximum:
    res.FBNS.FBNS(alg.it_tot,1) = max(res_aux(:));
end

if opc.flag_imposed == 2
    res.load.F_N(alg.it_tot,1)  = (sol.F_N - alg.load.F_N_aim)/ref.F_N;
end
end

function [fft2_Kernel,Kernel]  = construct_linear_Kernel(Nx1,Nx2,dx1,dx2,E_dash)
% Calculates the Kernel function for a linear convolution in the
% influence area of size Nx1*Nx2
% The Kernel is constructed for an imposed normal load pressure as explained by 
% Johnson, K. L., 2004. Contact mechanics. Cambridge: Cambridge University Press
% in equation (3.25) on P.54
% and the definition of 
% sld.E_dash = 2/((1 - sld.nu_low^2)/sld.E_low + (1 - sld.nu_up^2)/sld.E_up).
% The Kernel center is at the edge of the domain
% The Kernel domain is of size 2*Nx1*2*Nx2
% Input:
% Nx1                 [-]       Number of discretized points in x1-direction
% Nx2                 [-]       Number of discretized points in x2-direction
% dx1                 [m]       Length of discretized cell in x1-direction
% dx2                 [m]       Length of discretized cell in x2-direction
% E_dash              [Pa]      reduced modulus of elasticity
% Output:
% Kernel              [m/Pa]    Kernel
% fft2_Kernel         [?]       2-D fast Fourier transform of the Kernel
% -------------------------------------------------------------------------
Nx1_K = 2*Nx1;      % [-]       Number of discretized Kernel points in x1-direction
Nx2_K = 2*Nx2;      % [-]       Number of discretized Kernel points in x2-direction
dx1_mod = dx1/2;
dx2_mod = dx2/2;
% Determine distances:
% Code is written such that Nx1_K and Nx2_K can be odd or even:
% Indizes in x1-direction:
i_pos_end = floor(Nx1_K/2) + 1;
i_pos_sym = ceil(Nx1_K/2);
i_neg_sym = floor(Nx1_K/2) + 2;
i_pos     = 1:i_pos_end;
% Distances in x1-direction:
x1                  = zeros(1,Nx1_K);
x1(1:i_pos_end)     = dx1*(i_pos - 1);
x1(i_neg_sym:Nx1_K) = -flip(x1(2:i_pos_sym));
% Indizes in x2-direction:
j_pos_end = floor(Nx2_K/2) + 1;
j_pos_sym = ceil(Nx2_K/2);
j_neg_sym = floor(Nx2_K/2) + 2;
j_pos     = 1:j_pos_end;
% Distances in x2-direction:
x2                  = zeros(1,Nx2_K);
x2(1:j_pos_end)     = dx2*(j_pos - 1);
x2(j_neg_sym:Nx2_K) = -flip(x2(2:j_pos_sym));
% Convert vectors to matrizes:
[x1, x2] = ndgrid(x1,x2); 
% Construct Kernel:
term_1 = (x1 + dx1_mod).*log(ext_sqrt(x2 + dx2_mod, x1 + dx1_mod)./ext_sqrt(x2 - dx2_mod, x1 + dx1_mod));
term_2 = (x2 + dx2_mod).*log(ext_sqrt(x1 + dx1_mod, x2 + dx2_mod)./ext_sqrt(x1 - dx1_mod, x2 + dx2_mod));
term_3 = (x1 - dx1_mod).*log(ext_sqrt(x2 - dx2_mod, x1 - dx1_mod)./ext_sqrt(x2 + dx2_mod, x1 - dx1_mod));
term_4 = (x2 - dx2_mod).*log(ext_sqrt(x1 - dx1_mod, x2 - dx2_mod)./ext_sqrt(x1 + dx1_mod, x2 - dx2_mod));
clear x1;  clear x2; clear dx1_mod;  clear dx2_mod; 
Kernel        = 2/(pi*E_dash)*(term_1 + term_2 + term_3 + term_4);
clear term_1; clear term_2; clear term_3; clear term_4; 
fft2_Kernel   = fft2(Kernel);
function return_value = ext_sqrt(p,q)
% Auxiliary function for Kernel construction
return_value = p + sqrt(p.^2 + q.^2);
end
end

function [h_el] = compute_h_el(p,Nx1,Nx2,fft2_Kernel)
% Computes elastic displacement due to pressure field with linear
% convolution in Fourier space
% Input:
% p                 [Pa]        pressure field
% Nx1               [-]         number of discretized points in x1-direction
% Nx2               [-]       	number of discretized points in x2-direction
% fft2_Kernel       [?]         2-D fast Fourier transform of the Kernel
% Output:
% h_el              [m]         elastic deformation
% -------------------------------------------------------------------------
% Extend pressure field for linear convolution:
p_ext               = zeros(2*Nx1,2*Nx2);
p_ext(1:Nx1,1:Nx2)  = p(:,:);
% Compute convolution in Fourier space for better pefromance:
h_el_ext            = real(ifft2(fft2_Kernel.*fft2(p_ext)));
% Extract linear convolution:
h_el                = h_el_ext(1:Nx1,1:Nx2);
end

function [rho_l] = compute_rho_l(p_hd,fld,h)
if fld.rho_flag == 0
    % Constant denisty:
    rho_l = fld.rho_0*ones(h.Nx1,h.Nx2);
    
elseif fld.rho_flag == 1
    % Dowson-Higginson relation:
    % rho_l         [kg/m^3]    density of liquid lubricant
    % fld.rho_0     [kg/m^3]    liquid lubricant density at ambient pressure
    % p_hd          [Pa]        hydrodynamic pressure field
    rho_l = fld.rho_0.*(fld.dow_C1 + fld.dow_C2.*(p_hd - fld.p_cav))./(fld.dow_C1 + (p_hd - fld.p_cav));
    
end
end
function [prop] = compute_mu_l(sol,fld,prop,h)
if fld.mu_flag == 0
    % Constant dynamic viscosity:
    prop.mu_l = fld.mu_0*ones(h.Nx1,h.Nx2);
    
elseif fld.mu_flag == 1
    % Barus relation:
    % fld.barus_alpha     [1/Pa]      pressure viscosity coefficient, for mineral oils between 1e-8 and 2e-8 [1/Pa]
    prop.mu_l = fld.mu_0*exp(fld.barus_alpha*(sol.p_hd - fld.p_cav));
    
elseif fld.mu_flag == 2
    % Roelands relation:
    % mu_l          [Pas]       dynamic viscosity of liquid lubricant
    % fld.mu_0      [Pas]       liquid lubricant dynamic viscosity at ambient pressure
    % p_hd          [Pa]        hydrodynamic pressure field
    % fld.roe_z     [-]         pressure viscosity index
    % fld.roe_p0    [Pa]        constant in Roelands equation
    prop.mu_l = fld.mu_0.*exp((log(fld.mu_0) + 9.67)*(-1 + (1 + (sol.p_hd - fld.p_cav)/fld.roe_p0).^fld.roe_z));
    
end

% Consider shear thinning if fld.shear_thinning_model_flag > 0
if fld.shear_thinning_model_flag == 1
    % cut-off shear stress
    if fld.shear_thinning_limit_flag == 1
        % Consider constant limiting shear stress:
        % Use Couette shear stress as average limit along gap height for shear thinning:
        tau_av = abs(prop.u_r./h.h_m.*prop.mu_l);
        % Adjust effective dynamic viscosity:
        prop.mu_l(tau_av>fld.tau_max) = abs(fld.tau_max.*h.h_m(tau_av>fld.tau_max)/prop.u_r);
        
    elseif fld.shear_thinning_limit_flag == 2
        % Consider pressure dependent limiting shear stress:
        % Use Couette shear stress as average limit along gap height for shear thinning:
        tau_av = abs(prop.u_r./h.h_m.*prop.mu_l);
        % Adjust effective dynamic viscosity:
        prop.mu_l(tau_av>fld.tau_max_coeff*sol.p_hd) = abs(fld.tau_max_coeff*sol.p_hd(tau_av>fld.tau_max_coeff*sol.p_hd).*h.h_m(tau_av>fld.tau_max_coeff*sol.p_hd)/prop.u_r);
    
    end
elseif fld.shear_thinning_model_flag == 2
    % Eyring model
    gamma = abs(prop.u_r./h.h_m);                               % [1/s]     shear rate
    if fld.shear_thinning_limit_flag == 1
        % Consider constant limiting shear stress:
        % Adjust effective dynamic viscosity:
        prop.mu_l = fld.tau_max./gamma.*asinh(gamma.*prop.mu_l/fld.tau_max);
        
    elseif fld.shear_thinning_limit_flag == 2
        % Consider pressure dependent limiting shear stress:
        % Adjust effective dynamic viscosity:
        tau_max = fld.tau_max_coeff*sol.p_hd;
        prop.mu_l = tau_max./gamma.*asinh(gamma.*prop.mu_l./tau_max);
        
    end     
end
end

function [p_con,g,err] = elpl_contact_pressure_akchurin_linear(p_min,H,p_con_ini,z,h_s,Nx1,Nx2,fft2_Kernel,err_tol,it_max,h_ref)
% Erik Hansen, 26.08.2020
% Calculates the contact pressure occuring when a rigid smooth surface is
% loaded against an elastic half-space. The minimum and maximum contact pressures can be set to limiting values. Subsurface plastic flow within the half-space is not modelled.
% It is assumed that plastic flow only occures on the the half-space's surface if the maximum contact pressure is reached.
% The profile deformation is computed with the elastic half-space model using the Boundary Element Method (BEM).
% The code is very similar to and consists largely of the 
% algorithm described Akchurin et al., 2015:
% "On a model for the prediction of the friction coefficient in mixed lubrication based on a load-sharing concept with measured surface roughness"
% However, some changes were introduced in the treatment of non-contact and
% plastic contact points. Inspiration was taken from Vollebregt, 2014:
% "A new solver for the elastic normal contact problem using conjugate gradients, deflation, and an FFT-based preconditioner"
% The Matlab implementation of the linear convolution in the fourier space
% was largely inspired by Sainsot and Lubrecht, 2011:
% "Efficient solution of the dry contact of rough surfaces: a comparison of fast Fourier transform and multigrid methods"
% Input:
% p_min             [Pa]    pressure in non contact zones as a lower limit for the contact pressure
% H                 [Pa]    hardness of the softer material as a maximum limit of the contact pressure
% p_con_ini         [Pa]    initial guess of contact pressure field
% z                 [m]     undeformed profile
% h_s               [m]     seperating distance between upper and lower reference planes
% Nx1               [-]     number of points in x1-direction
% Nx2               [-]     number pf points in x2-direction
% fft2_Kernel       [?]     2-D fast Fourier transform of the Kernel
% err_tol           [-]     relative error tolerance
% it_max            [-]     maximum number of iterations
% h_ref             [m]     reference length of relative error
% Output:
% p_con             [Pa]    contact pressure field
% g                 [m]     residual of the gap height distribution
% err               [-]     relative error
% -------------------------------------------------------------------------
G_ref  = h_ref^2;                               % [m^2] reference norm of the residual
p_con = p_con_ini;                              % [Pa] pressure field
[u] = compute_h_el(p_con,Nx1,Nx2,fft2_Kernel);  % [m] elastic deformation 
g = -u + z - h_s;                               % [m] residual of the gap height distribution

% Find indices of points that are in the non-contact, elastic and plastic
% domain due to the pressure distribution. Non-contact and plastic points 
% are also evaluated whether they are correctly or not in surface contact
% due to the gap height distribution
A_el    = find(p_con>p_min  &p_con<H);
A_nc_cr = find(p_con<=p_min &g<0);
A_nc_wr = find(p_con<=p_min &g>=0);
A_pl_cr = find(p_con>=H     &g>=0);
A_pl_wr = find(p_con>=H     &g<0);
% Within these points, the pressure distribution needs to be adjusted:
A_free = union(union(A_el,A_nc_wr),A_pl_wr);
 
G = sum(g(A_free).*g(A_free));  % [m^2] norm of the residual
G_old = 1;                      % [m^2] previous norm of the residual
delta = 0;                      % [-] flag whether to use conjugate gradient or steepest descend
err = zeros(it_max,1);          % [-] relative error
i_it=0;                         % [-] iteration counter        
t = zeros(Nx1,Nx2);             % [Pa] search direction
while i_it == 0 || (err(i_it)>err_tol && i_it<=it_max)
    i_it = i_it + 1;    
    
    % Find search direction:
    t(A_free)   = g(A_free) + delta*(G/G_old)*t(A_free);
    t(A_nc_cr)  = 0;
    t(A_pl_cr)  = 0;
    clear delta;
    
    % Determine step length tau:
    [t_conv]    = compute_h_el(t,Nx1,Nx2,fft2_Kernel);     
    r           = -t_conv;
    clear t_conv;
    tau = (sum(g(A_free).*t(A_free)))/(sum(r(A_free).*t(A_free)));
    clear r;

    % Update pressure:
    p_con(A_free)       = p_con(A_free) - tau*t(A_free);
    p_con(p_con<p_min)  = p_min;
    p_con(p_con>H)      = H;
    
    % Compute elsatic deformation to find new residual of the gap height distribution
    [u] = compute_h_el(p_con,Nx1,Nx2,fft2_Kernel); 
    g = -u + z - h_s;
    clear u;
    
    % Find indices of points that are in the non-contact, elastic and plastic domain:
    A_el    = find(p_con>p_min  &p_con<H);
    A_nc_cr = find(p_con<=p_min &g<0);
    A_nc_wr = find(p_con<=p_min &g>=0);
    A_pl_cr = find(p_con>=H     &g>=0);
    A_pl_wr = find(p_con>=H     &g<0);
    % Within these points, the pressure distribution needs to be adjusted:
    A_free = union(union(A_el,A_nc_wr),A_pl_wr);
    
    % Determine whether to use conjugate gradient or steepest descend in the next ieration
    if isempty(A_nc_wr) && isempty(A_pl_wr)
    	delta = 1;
    else
        delta = 0;
    end

    % Save G for the next iteration:
    G_old = G;
    % Compute norm of the residual:
    G = sum(g(A_free).*g(A_free));
    % Compute relative error
    err(i_it) = abs(G - G_old)/G_ref;
end
% Resize err:
err = err(1:i_it,1);
end

function [tau_hd_low,tau_hd_up,dp_hd_dx1] = compute_tau_hd(h,sol,prop)
% This function computes the shear stresses on the upper and lower surface
% and the pressure gradient:
dp_hd_dx1               = zeros(h.Nx1,h.Nx2);
i                       = 2:h.Nx1-1;
j                       = 1:h.Nx2;
dp_hd_dx1(i,j)          = (sol.p_hd(i+1,j)      - sol.p_hd(i-1,j)    )/(2*h.dx1);               % [Pa/m]    hydrodynamic pressure gradient in x1-direction
dp_hd_dx1(1,j)          = (sol.p_hd(2,j)        - sol.p_hd(1,j)      )/h.dx1;
dp_hd_dx1(h.Nx1,j)      = (sol.p_hd(h.Nx1,j)    - sol.p_hd(h.Nx1-1,j))/h.dx1;
% Assuming mu does not depend on the cavity fraction: mu = mu_l:
tau_hd_up               = ( 1/2*h.h_m.*dp_hd_dx1 + prop.mu_l.*(1 - sol.thet).*prop.u_r./h.h_m);     % [Pa]  hydrodynamic shear stress \tau_{31} at the upper surface 
tau_hd_low              = (-1/2*h.h_m.*dp_hd_dx1 + prop.mu_l.*(1 - sol.thet).*prop.u_r./h.h_m);     % [Pa]  hydrodynamic shear stress \tau_{31} at the lower surface 
end
