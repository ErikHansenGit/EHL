close all; clc; clearvars -except autorun;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EHL setup code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main developer:
% Erik Hansen: erik.hansen@kit.edu
% Further developers who contributed to the code:
% Altay Kacan, 08.01.2021: Implementation of the unsteady terms in the
% Reynolds equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information:
% This code creates and saves a setup to be then used by the script
% EHL_02_mainprocess.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation Identifier:
alg.Sim_ID = autorun.alg.Sim_ID{autorun.sim.it};

%% Fluid settings:
fld.h_min       = 1e-9;  
fld.p_amb       = 1e5;                      % [Pa]      ambient pressure, used for load force computation
fld.p_cav       = 0e0;                      % [Pa]      cavitation pressure

% Boundary conditions:
% Hydrodynamic pressure:
fld.flag_p_bound    = 0;                    % [-]       flag, which boundary conditions for the hydrodynamic pressure are used: 
% 0) Dirichlet at West, North, South and East boundary
% 1) Dirichlet at West and East boundary, Neumann at North and South boundary
% 2) Dirichlet at West boundary, Neumann at North, South and East boundary
if fld.flag_p_bound == 0
    fld.p_S_side      = fld.p_amb;          % [Pa]  pressure on the South side, must be in interval of [fld.p_cav,inf]
    fld.p_W_side      = fld.p_amb;          % [Pa]  pressure at the inlet, must be in interval of [fld.p_cav,inf] 
    fld.p_E_side      = fld.p_amb;          % [Pa]  pressure at the outlet, must be in interval of [fld.p_cav,inf]  
    fld.p_N_side      = fld.p_amb;          % [Pa]  pressure on the North side, must be in interval of [fld.p_cav,inf]  
elseif fld.flag_p_bound == 1
    fld.p_S_side      = NaN;                % [Pa]  pressure on the South side, must be in interval of [fld.p_cav,inf]
    fld.p_W_side      = fld.p_amb;          % [Pa]  pressure at the inlet, must be in interval of [fld.p_cav,inf] 
    fld.p_E_side      = fld.p_amb;          % [Pa]  pressure at the outlet, must be in interval of [fld.p_cav,inf]  
    fld.p_N_side      = NaN;                % [Pa]  pressure on the North side, must be in interval of [fld.p_cav,inf]
elseif fld.flag_p_bound == 2
    fld.p_S_side      = NaN;                % [Pa]  pressure on the South side, must be in interval of [fld.p_cav,inf]
    fld.p_W_side      = fld.p_amb;          % [Pa]  pressure at the inlet, must be in interval of [fld.p_cav,inf] 
    fld.p_E_side      = NaN;                % [Pa]  pressure at the outlet, must be in interval of [fld.p_cav,inf]  
    fld.p_N_side      = NaN;                % [Pa]  pressure on the North side, must be in interval of [fld.p_cav,inf]   
end
    
% Cavity fraction:
fld.flag_thet_bound    = 0;                 % [-]       flag, which boundary conditions for the hydrodynamic pressure are used: 
% 0) Dirichlet at West, North, South and East boundary
% 1) Dirichlet at West and East boundary, Neumann at North and South boundary
% 2) Dirichlet at West boundary, Neumann at North, South and East boundary
if fld.flag_thet_bound == 0
    fld.thet_S_side      = 0.0;          % [-]  cavity fraction on the South side, must be in interval of [0,1]
    fld.thet_W_side      = 0.0;          % [-]  cavity fraction at the inlet, must be in interval of [0,1] 
    fld.thet_E_side      = 0.0;          % [-]  cavity fraction at the outlet, must be in interval of [0,1]  
    fld.thet_N_side      = 0.0;          % [-]  cavity fraction on the North side, must be in interval of [0,1]  
elseif fld.flag_thet_bound == 1
    fld.thet_S_side      = NaN;          % [-]  cavity fraction on the South side, must be in interval of [0,1]
    fld.thet_W_side      = 0.0;          % [-]  cavity fraction at the inlet, must be in interval of [0,1] 
    fld.thet_E_side      = 0.0;          % [-]  cavity fraction at the outlet, must be in interval of [0,1]  
    fld.thet_N_side      = NaN;          % [-]  cavity fraction on the North side, must be in interval of [0,1]  
elseif fld.flag_thet_bound == 2
    fld.thet_S_side      = NaN;          % [-]  cavity fraction on the South side, must be in interval of [0,1]
    fld.thet_W_side      = 0.0;          % [-]  cavity fraction at the inlet, must be in interval of [0,1] 
    fld.thet_E_side      = NaN;          % [-]  cavity fraction at the outlet, must be in interval of [0,1]  
    fld.thet_N_side      = NaN;          % [-]  cavity fraction on the North side, must be in interval of [0,1]     
end


% Density:
fld.rho_0       = 850;                                                      % [kg/m^3]  lubricant density at ambient pressure, value does not actually matter since dimensionless Reynolds equation is solved
fld.rho_flag    = 1;                                                        % [-]       flag whether liquid lubricant density is constant (0) or Dowson-Higginson relation is used for pressure dependency (1)
if fld.rho_flag == 1
    fld.dow_C1  = 2.22e9;                                                   % [Pa]      constant in Dowson Higginson relation
    fld.dow_C2  = 1.66;                                                     % [-]      coefficient in Dowson Higginson relation
end


% Dynamic viscosity:
fld.mu_0        = 1.0e-2;                                                   % [Pas]     lubricant dynamic viscosity at ambient pressure
fld.mu_flag     = 1;                                                        % [-]       flag whether liquid lubricant dynamic viscosity is constant (0), Barus relation (1) or Roelands relation (2) is used for pressure dependency 
if fld.mu_flag == 1
    fld.barus_alpha = 1.2e-8;                                               % [1/Pa]    pressure viscosity coefficient, for mineral oils between 1e-8 and 2e-8 [1/Pa]
elseif fld.mu_flag == 2
    fld.roe_alpha   = NaN;                                                  % [1/Pa]    pressure viscosity coefficient, for mineral oils between 1e-8 and 2e-8 [1/Pa]
    fld.roe_p0      = NaN;                                                  % [Pa]      constant in Roelands equation
    fld.roe_z       = fld.roe_alpha*fld.roe_p0/(log(fld.mu_0) + 9.67);      % [-]       pressure viscosity index
end


% Shear thinning:
fld.shear_thinning_model_flag = 0;                                          % [-]       flag how shear thinning is considered
% 0: no shear thinning
% 1: cut-off
% 2: Eyring
if fld.shear_thinning_model_flag > 0    
    fld.shear_thinning_limit_flag = 1;                                      % [-]       flag how shear thinning limit is considered
    if fld.shear_thinning_limit_flag == 1
        % Constant limit
        fld.tau_max     = 5e6;                                              % [Pa]      limiting shear stress
    elseif  fld.shear_thinning_limit_flag == 2
        % Varying limit
        fld.tau_max_coeff = NaN;                                            % [-]       fraction of hydrodynamic pressure that is used for the limiting shear stress
    end
end

%% Solid settings:
sld.flag_elastic    = autorun.sld.flag_elastic(autorun.sim.it);             % [-]   flag whether solid body is elastic on macroscopic scale (true) or rigid (false)
if sld.flag_elastic
    sld.nu_low      = 0.3;                                                  % [-]   Poisson ratio of lower body
    sld.E_low       = 210e9;                                                % [Pa]  Young's modulus of lower body
    sld.nu_up       = 0.3;                                                  % [-]   Poisson ratio of upper body
    sld.E_up        = 210e9;                                                % [Pa]  Young's modulus of upper body
    sld.E_dash  = 2/((1 - sld.nu_low^2)/sld.E_low + (1 - sld.nu_up^2)/sld.E_up);  % [Pa]   reduced modulus of elasticity
end

%% Operating condition:
opc.flag_imposed      = 1;                                                  % [-] flag whether a constant rigid body displacement h_d (1) or a constant normal load force W (2) is imposed
if opc.flag_imposed == 1
    opc.h_d           = 0e-0;                                               % [m] imposed rigid body displacement between upper and lower surface
elseif opc.flag_imposed == 2
    opc.F_N           = NaN;                                                % [N] imposed normal load force in negative x3-direction     
end
% To specify several operating conditions, just set up opc.u_up_ini and/or
% opc.u_low_ini as 1D arrays, e.g.: opc.u_up_ini = [9.0e-2 1.0e-2];                                          
opc.u_up_ini            = 1;                                % [m/s] velocity of upper body in x1-direction, the disk
opc.u_low_ini           = 0;                                % [m/s] velocity of lower body in x1-direction, the ball

[opc.u_up, opc.u_low]   = ndgrid(opc.u_up_ini,opc.u_low_ini);             % constructing matrices that correspond to each combination of the given initial value

opc.u_r                 =  opc.u_up - opc.u_low;                          % [m/s] relative velocities for each operating condition
opc.u_m                 = (opc.u_up + opc.u_low)/2;                       % [m/s] mean velocities for each operating condition, must be larger than 0 
opc.srr                 = (opc.u_low - opc.u_up)/(opc.u_low + opc.u_up);  % [-] slide-to-roll ratio // -1 for u_low=0 // 0 for u_low=u_up 

opc.N                   = numel(opc.u_m);                                 % [-] number of of operating conditions

%% Geometry settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Study  B: Single Dimple - Bertocchi et al., 2013:
% Input:
% Gap height:
geo.h_min           = 1.0e-6;       % [m]       minimum gap height
geo.h_max           = 1.1e-6;       % [m]       maximum gap height
geo.h_p             = 0.4e-6;       % [m]       dimple depth

% x1-direction:
geo.L_x1            = 2e-2;                             % [m]       domain lenght
geo.l_x1            = 6.0e-3;                           % [m]       dimple length
geo.l_p_x1          = 4.0e-3;                           % [m]       distance between dimple and inlet
geo.Nx1             = 2^7;                              % [-]       number of discretization cells
geo.x1              = linspace(0,geo.L_x1,geo.Nx1);     % [m]       coordinate of cell centers
geo.dx1             = geo.L_x1/(geo.Nx1 - 1);           % [m]       spacing

% x2-direction:
geo.L_x2            = autorun.geo.L_x2(autorun.sim.it); % [m]       domain width
geo.l_x2            = autorun.geo.l_x2(autorun.sim.it); % [m]       dimple width
geo.l_p_x2          = (geo.L_x2 - geo.l_x2)/2;          % [m]       distance between dimple and boundary
geo.Nx2             = 2^6+1;                            % [-]       number of discretization cells
geo.x2              = linspace(0,geo.L_x2,geo.Nx2);     % [m]       coordinate of cell centers
geo.dx2             = geo.L_x2/(geo.Nx2 - 1);           % [m]       spacing

% Computation:
[x1_matr,x2_matr] = ndgrid(geo.x1,geo.x2);
geo.h_g_ma  = geo.h_max + (geo.h_min - geo.h_max)/geo.L_x1*x1_matr;

cond_x1_ins = x1_matr>geo.x1(1)+geo.l_p_x1 & x1_matr<geo.x1(1)+geo.l_p_x1+geo.l_x1; 
cond_x2_ins = x2_matr>geo.x2(1)+geo.l_p_x2 & x2_matr<geo.x2(1)+geo.l_p_x2+geo.l_x2;
clear x1_matr; clear x2_matr;

cond_ins = cond_x1_ins==1 & cond_x2_ins==1;
clear cond_x1_ins; clear cond_x2_ins;

geo.h_g_ma(cond_ins) = geo.h_g_ma(cond_ins) + geo.h_p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save profile information:
geo.x3     = -geo.h_g_ma;                                          % [m]   macroscopic profile
geo.N      = geo.Nx1*geo.Nx2;

%% Unsteady geometry variation:
% Transient algorithm settings:
alg.flag_unsteady = false;               % [-]   boolean, whether simulation is unsteady (true) or steady (false) 

if alg.flag_unsteady
    %Define micro-texture parameters:
    geo_Ti.x1_0         =  geo.x1(1);             % [m] initial x1 coordinate of the micro-texture center
    geo_Ti.x2_0         =  0;                     % [m] initial x2 coordinate of the micro-texture center
    geo_Ti.u_tex        =  opc.u_low;             % [m/s] velocity of the micro-texture 

    delta_t_nd  = NaN;                                                % [-]  dimensionless time step size, same size as dimensionless space resolution same as Mourier (2006)
    opc.delta_t = NaN;                                 % [s]  time step size for each operating condition 
    opc.N_t     = round((geo.x1(geo.Nx1) - geo.x1(1))./(opc.delta_t.*geo_Ti.u_tex) + 1);    % [-]   amount of discrete time steps 
    
    % Initialize the gap height change and the matrix containing the dimple center coordinates:
    time        = zeros(max(opc.N_t(:)),opc.N);                  % [s] discrete times
    h_t_ma      = zeros(geo.Nx1,geo.Nx2,max(opc.N_t(:)),opc.N);  % [m] 4 dimensional array (2 spatial, 1 temporal, 1 extra dimension for operating condition)
    x1_c_vec    = zeros(max(opc.N_t(:)),opc.N);                  % [m] 2 dimensional array x1 coordinates of the dimple center for each time point and operating condition 
    x2_c_vec    = zeros(max(opc.N_t(:)),opc.N);                  % [m] 2 dimensional array x2 coordinates of the dimple center for each time point and operating condition 
    
    % Loop for going through all operating conditions:
    for it_OC = 1:opc.N 

        % Time loop to precalculate imposed transient gap height change for all time points:
        for it_time = 1:opc.N_t(it_OC)

             fprintf("Precalculating time point: %i/%i for operating condition: %i/%i\n",it_time,opc.N_t(it_OC),it_OC,opc.N);
             % Calculate corresponding time:
             time(it_time,it_OC)  = (it_time-1)*opc.delta_t(it_OC);     % At the first time point time is = 0, the time step size depends on the operating condition  
             % Calculate gap height change due to micro-texture for this time point:
             x1_c        = geo_Ti.x1_0 + geo_Ti.u_tex(it_OC)*time(it_time,it_OC);      % [m] time dependent x1 coordinate of the texture center    
             x1_c_vec(it_time,it_OC) = x1_c;                            % [m] saved value of the x1 coordinate of the dimple center  
             x2_c        = 0;                                           % [m] x2 coordinate of the texture center
             x2_c_vec(it_time,it_OC) = x2_c;                            % [m] saved value of the x1 coordinate of the dimple center  
             h_t_ma(:,:,it_time,it_OC)   = NaN;                         % [m] gap height influence due to the micro-texture, the third dimension is time and the fourth dimension is the corresponding operating condition
        end  
    end
    
else
    opc.N_t = ones(size(opc.u_m,1),size(opc.u_m,2));% [-]   amount of discrete time steps 
end

%% Algorithm settings:
% Genereal relative tolerance:
alg.it_tot_max          = 1e3;          % [-]   maximum number of total iterations until simulation is aborted
alg.FBNS.toll           = 1e-6;         % [-]   tolerance for FBNS residual checks
alg.FBNS.alpha_p_nd     = 5e-2;         % [-]   relaxation factor for update of hydrodynmic pressure
alg.FBNS.alpha_thet     = 5e-2;         % [-]   relaxation factor for update of cavity fraction
alg.FBNS.flag_enforce_limit_p_nd_red_min = true;  % [-]   boolean whether to truncate updated p_nd_red_min at 0 (true) or not (false)
alg.FBNS.flag_enforce_limit_thet_red_min = true;  % [-]   boolean whether to truncate updated thet_red_min at 0 (true) or not (false)
alg.FBNS.flag_enforce_limit_p_red_max    = false; % [-]   boolean whether to truncate updated p_nd_red_min at alg.FBNS.limit_p_nd_red_max (true) or not (false)
if alg.FBNS.flag_enforce_limit_p_red_max 
    alg.FBNS.limit_p_red_max = NaN;               % [Pa]   maximum limit for the hydrdynamic pressure    
end
alg.FBNS.flag_enforce_limit_thet_red_max    = true;  % [-]   boolean whether to truncate updated thet_red_min at 1 (true) or not (false)

alg.FBNS.flag_res       = 10;           % [-]   flag, which residual definition is used 
% The FBNS residual defintions can be found and added in the function FBNS_det_residual
% When in doubt choose alg.FBNS.flag_res = 10 (strict with respect to, G,F,p_hd and theta)
%                   or alg.FBNS.flag_res = 11 (strict with respect to, G,p_hd)
alg.FBNS.flag_print_it  = 1;            % [-]   flag, whether information about each iteration of the FBNS algorithm is printed (1) or not (0)

alg.FBNS.flag_Couette_Discrete = autorun.alg.FBNS.flag_Couette_Discrete(autorun.sim.it);     % [-]   flag, which discretization scheme is used for the Couette term 
% (maximum order is 2 since midpoint rule is always used for integral flux approximations): 
% 1) 1.Order Upwind, 2) 2.Order QUICK, 3) 2.Order LUI, 4) 2.Order CUI

alg.FBNS.flag_A_h_S  = true;     % [-]   boolean whether (true) or not (false) to use the South Diagonal of A_h: true always advised
alg.FBNS.flag_A_h_WW = false;    % [-]   boolean whether (true) or not (false) to use the WestWest Diagonal of A_h: true advised if LUI is used, otherwise false is advised
alg.FBNS.flag_A_h_E  = true;     % [-]   boolean whether (true) or not (false) to use the East Diagonal of A_h: true always advised
alg.FBNS.flag_A_h_N  = true;     % [-]   boolean whether (true) or not (false) to use the North Diagonal of A_h: true always advised

alg.FBNS.flag_p_hd_ini = 1;      % [-]   flag, which initial guess is used for the hydrodynamic pressure distribution
% 1) fld.p_amb, 2) elastic dry contact solution
if alg.FBNS.flag_p_hd_ini == 2
    alg.FBNS.p_hd_ini.toll      = NaN;              % [-]   tolerance of during computation of inital pressure distribution
    alg.FBNS.p_hd_ini.it_max    = NaN;              % [-]   maximum number of allowed iteration steps before computation of inital pressure distribution is stopped
end

% Reference quantities:
alg.FBNS.ref.x1     = geo.L_x1;                     % [m]
alg.FBNS.ref.x2     = geo.L_x2;                     % [m]
alg.FBNS.ref.h      = geo.h_min;                    % [m]
alg.FBNS.ref.mu     = fld.mu_0;                     % [Pas]
alg.FBNS.ref.rho    = fld.rho_0;                    % [kg/m^3]
alg.FBNS.ref.u      = opc.u_m;                      % [m/s]
alg.FBNS.ref.p      = 10e6;                         % [Pa]
if alg.flag_unsteady
    alg.FBNS.ref.t  = NaN;                          % [s]
end

alg.load.toll = alg.FBNS.toll;                      % [-]   tolerance of load balance eqaution
if opc.flag_imposed == 2
    
    alg.load.ref.F_N    = NaN;                      % [N]
    
    alg.load.h_d_ma_ini = NaN;                      % [m]   initial guess for rigid body displacement

    alg.load.K_P = NaN;                             % [-]   Propotional coefficient of load balance PID-controller
    alg.load.K_I = NaN;                             % [-]   Integration coefficient of load balance PID-controller
    alg.load.K_D = NaN;                             % [-]   Derivative  coefficient of load balance PID-controller
end

% Advanced algorithm settings:
alg.flag_profile_viewer = false;            % [-]   boolean, whether profile viewer is activated (true) or not (false)
alg.flag_limit_num_comp_threads  = false;   % [-]   boolean, whether whther the number of computational threads is limited (true) or not (false)
if alg.flag_limit_num_comp_threads
    alg.max_num_comp_threads = 8;           % [-]   maximum number of computational threads
end

%% Save setup data:
% Define path:
output_path = './../data/EHL_02_mainprocess/Input/';
% Create folder/ delete previous contents:
warning('off','MATLAB:MKDIR:DirectoryExists')
mkdir (output_path)
rmdir(output_path, 's')        % remove old setups if they exist
mkdir (output_path)
% Save data:
save (fullfile(output_path,'fld.mat'),'fld');
save (fullfile(output_path,'sld.mat'),'sld');
save (fullfile(output_path,'geo.mat'),'geo');
save (fullfile(output_path,'opc.mat'),'opc');
save (fullfile(output_path,'alg.mat'),'alg');

if alg.flag_unsteady
    save (fullfile(output_path,'geo_Ti.mat'),'geo_Ti');
    for it_OC=1:opc.N 
        % Create directories for each operating condition:
        output_path_OC   = fullfile(output_path,sprintf('h_time/OC_%i/',it_OC));
        mkdir (output_path_OC)
        for it_time=1:opc.N_t(it_OC)
            % Create directories for each time point:    
            output_path_time = fullfile(output_path_OC,sprintf('Time_%i/',it_time));
            mkdir (output_path_time)
            % Create and save the corresponding slices of the 4D array:
            slc.h_t_ma = h_t_ma(:,:,it_time,it_OC);
            slc.x1_c   = x1_c_vec(it_time,it_OC);
            slc.x2_c   = x2_c_vec(it_time,it_OC);
            slc.t      = time(it_time,it_OC);
            save(fullfile(output_path_time,'slc.mat'),'slc');
        end
    end
end

%% Plot settings:
flag_plot_setup_data = false;
if flag_plot_setup_data 
    
    % Line plots:
    KIT_colorlist={[0,150,130]/255,[162 34 35]/255,[70 100 170]/255,[252 229 0]/255,[140 182 60]/256,[223 155 27]/255,[167 130 46]/255,[163 16 124]/255,[35 161 224]/255};
    % Surface plots:
    colmap = parula;
    % Size:
    widthlines          = 1.2;
    sizeoffonts         = 11;
    sizeoflegendfonts   = 11;

    % Plots:
    % Surface plot gap height:
    fig = figure('Units','centimeters','Position',[10 08 7 7]);
    if alg.flag_unsteady
        surf(geo.x1,geo.x2,geo.h_g_ma'+h_t_ma(:,:,1,1)')
    else
        surf(geo.x1,geo.x2,geo.h_g_ma')
    end
    xlabel('x_1 [m]','fontsize',sizeoffonts);
    ylabel('x_2 [m]','fontsize',sizeoffonts);
    zlabel('h_m [m]','fontsize',sizeoffonts);
    title('Mean gap height','fontsize',sizeoffonts)
    fig.CurrentAxes.ZDir = 'Reverse';
    shading interp
    colormap(flipud(colmap))
    light
    lighting 'gouraud'
    lightangle(-90,-30)
    clear fig;
end

fprintf('\n--------------------------------');
fprintf('\n<strong>Setup finished!</strong>');
fprintf('\nSimulation ID: ''%s''',alg.Sim_ID);
fprintf('\n--------------------------------\n');
