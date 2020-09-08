close all; clc; clear all;
% Setup of the EHL solver settings
%
% This script exports the data specified in this script in a format 
% which can be used by the script
% "EHL_02_mainprocess.m" to solve the EHL problem.
% The exact file paths of the output data can be specified in 
% the "File path information and output" section of this script.
%
% This code replicates the setup of:
% Křupka, I., Hartl, M., Urbanec, L., & Čermák, J. (2007). 
% Single dent within elastohydrodynamic contact-comparison between experimental and numerical results. 
% Proceedings of the Institution of Mechanical Engineers, Part J: Journal of Engineering Tribology, 221(6), 635-644.
% 
% Erik Hansen, 07.09.2020

% Set input information:
%% Fluid settings:
fld.p_amb       = 101325;                                                   % [Pa]      ambient pressure
fld.p_cav       = fld.p_amb;                                                % [Pa]      cavitation pressure
fld.rho_0       = 850;                                                      % [kg/m^3]  lubricant density at ambient pressure
fld.mu_0        = 0.320;                                                    % [Pas]     lubricant dynamic viscosity at ambient pressure
fld.alpha       = 31e-9;                                                    % [1/Pa]    pressure viscosity coefficient 
fld.tau_max     = 1e7;                                                      % [Pa]      limiting shear stress
fld.h_min       = 1e-9;                                                     % [m]       minimum macroscopic gap height to prevent gap height from becoming zero

%% Solid seetings:
sld.nu_low      = 0.3;                                                      % [-]   Poisson ratio of lower body
sld.E_low       = 210e9;                                                    % [Pa]  Young's modulus of lower body
sld.nu_up       = 0.208;                                                    % [-]   Poisson ratio of upper body
sld.E_up        = 81e9;                                                     % [Pa]  Young's modulus of upper body
sld.E_dash      = 2/...
    ((1 - sld.nu_low^2)/sld.E_low + (1 - sld.nu_up^2)/sld.E_up);            % [Pa]   reduced modulus of elasticity

%% Operating condition:
opc.W           = 29;                                                       % [N] Normal load force in negative x3-direction
opc.u_low       = 0;                                                        % [m/s] velocity of lower body in x1-direction
opc.u_up        = [2*0.0342];                                               % [m/s] velocity of upper body in x1-direction
opc.N           = length(opc.u_up);                                         % [-]   number of of operating conditions

%% Geometry settings:
% Define geometry parameters
geo.Nx1         = 2^8;                                                      % [-]   number of discretization points in the x1-direction
geo.Nx2         = 2^8;                                                      % [-]   number of discretization points in the x2-direction
x1_Wb           = -0.42e-3;                                                 % [m]   x1-coordinate of cell at the West boundary
x1_Eb           = 0.25e-3;                                                  % [m]   x1-coordinate of cell at the East boundary
x2_Sb           = -0.25e-3;                                                 % [m]   x2-coordinate of cell at the South boundary
x2_Nb           = 0.25e-3;                                                  % [m]   x2-coordinate of cell at the North boundary
geo.Rx1         = 1.27e-2;                                                  % [m]   ball radius

% Construct geometry:
geo.dx1         = (x1_Eb - x1_Wb)/(geo.Nx1 - 1);                            % [m]   spacing in x1-direction
geo.dx2         = (x2_Nb - x2_Sb)/(geo.Nx2 - 1);                            % [m]   spacing in x2-direction
geo.x1          = linspace(x1_Wb,x1_Eb,geo.Nx1);                            % [m]   x1-coordinates with uniform discretization
geo.x2          = linspace(x2_Sb,x2_Nb,geo.Nx2);                            % [m]   x2-coordinates with uniform discretization
[x1_matr, x2_matr] = ndgrid(geo.x1,geo.x2);
% Ball-on-disc tribometer:
geo.h_prof_ma   = (x1_matr.^2)/(2*geo.Rx1) + (x2_matr.^2)/(2*geo.Rx1);      % [m]   gap height variation induced by macroscopic profile
geo.prof_ma     = -geo.h_prof_ma;                                           % [m]   macroscopic profile
clear x1_matr; clear x2_matr;

%% Algorithm seetings:
% Genereal relative tolerance:
alg.toll                = 1e-5;         % [-]   relative tolerance for residual checks
% Secant algorithm to solve load balance equation:
alg.it_max_W            = 3e2;          % [-]   maximum number of iterations for checking of load balance equation
alg.h_d_ma_nm2          = -2.1e-6;      % [m]   first initial guess of rigid body displacement in secant algorithm
alg.h_d_ma_nm1          = -2.0e-6;      % [m]   second initial guess of rigid body displacement in secant algorithm
alg.h_d_ma_max          = -1.5e-6;      % [m]   upper limit of rigid body displacement in secant algorithm
alg.h_d_ma_min          = -2.1e-6;      % [m]   lower limit of rigid body displacement in secant algorithm
alg.h_delta             = 1e-7;         % [m]   step size for rigid body displacement when limit is reached in secant algorithm to reset the second initial guess of rigid body displacement
alg.toll_W_delta_p      = 1e-1;         % [-]   relative FBNS pressure residual which has to be reached in order to evaluate the load balance equation
% Contact pressure algorithm to find an initial guess for the hydrodynamic pressure:
alg.it_max_con          = 1e2;          % [-]   maximum number of iterations in contact pressure algorithm
alg.h_ref_con           = 1e-8;         % [m]   reference length of relative residual in contact pressure algorithm
% FBNS algorithm to solve the EHL problem for hydrodynamic pressure and cavity fraction:
alg.it_max_FBNS         = 2e1;          % [-]   maximum number of iterations until load balance equation is checked
alg.p_ref_FBNS          = opc.W/...
    (geo.Nx1*geo.dx1*geo.Nx2*geo.dx2);  % [Pa]  reference pressure for relative FBNS pressure residual
alg.alpha_FBNS_above    = 5e-2;         % [-]   relaxation factor for FBNS algorithm, used in first FBNS iteration or when the relative FBSN residual is greater than alg.res_alpha_FBNS    
alg.alpha_FBNS_below    = 5e-2;         % [-]   relaxation factor for FBNS algorithm, when the relative FBSN residual is smaller than alg.res_alpha_FBNS   
alg.res_alpha_FBNS      = 1e-1;         % [-]   boundary for the relative FBNS residual to switch between relaxation factors   

% Performace settings:
alg.flag_profile_viewer = false;        % [-]   enable profiler and display results in the end
alg.max_num_com_threads = 8;            % [-]   maximum number of computational threads

%% File path information and output:
% Output:
% Output path of setup information(->Input path of EHL_02_mainprocess.m):
output_path     = sprintf('%s','./../data/EHL_02_mainprocess/Input');
mkdir (output_path)
% Output path of solver results(->Output path of EHL_02_mainprocess.m):
inf.output_path_mainprocess     = sprintf('%s','./../data/Krupka/EHL_02_mainprocess/Output');
mkdir (output_path)
% Save output:
save(fullfile(output_path,'/alg.mat'),'alg');
save(fullfile(output_path,'/fld.mat'),'fld');
save(fullfile(output_path,'/sld.mat'),'sld');
save(fullfile(output_path,'/geo.mat'),'geo');
save(fullfile(output_path,'/opc.mat'),'opc');
save(fullfile(output_path,'/inf.mat'),'inf');

%% Plot:
% Plot settings:
% Line plots:
KIT_colorlist={[0,150,130]/255,[162 34 35]/255,[70 100 170]/255,[252 229 0]/255,[140 182 60]/256,[223 155 27]/255,[167 130 46]/255,[163 16 124]/255,[35 161 224]/255};
% Surface plots:
colmap = parula;
% Size:
widthlines          = 1.2;
sizeoffonts         = 11;
sizeoflegendfonts   = 11;

% Macroscopic geometry:
figure('Units','centimeters','Position',[00 00 10 10])
surf(geo.x1,geo.x2,geo.prof_ma')
xlabel('x_1 [m]','fontsize',sizeoffonts);
ylabel('x_2 [m]','fontsize',sizeoffonts);
zlabel('x_3 [m]','fontsize',sizeoffonts);
title('Pin profile','fontsize',sizeoffonts)
shading interp
material dull
colormap(colmap)
camlight

% Macroscopic gap height:
[x2_0,x2_0_i] = min(abs(geo.x2));
fig = figure('Units','centimeters','Position',[11 00 10 10]);
plot(geo.x1*1e6,geo.h_prof_ma(:,x2_0_i)*1e9,'LineWidth',widthlines,'Color',KIT_colorlist{1})
grid on
xlabel('x_1 [\mu m]','fontsize',sizeoffonts);
ylabel('h_{ma,ad} [nm]','fontsize',sizeoffonts);
xlim([geo.x1(1)*1e6 geo.x1(geo.Nx1)*1e6])
title("x_2= " + x2_0*1e6 + "\mu m",'fontsize',sizeoffonts)
clear x2_0; clear x2_0_i;

fprintf('\nSetup successful! \n')