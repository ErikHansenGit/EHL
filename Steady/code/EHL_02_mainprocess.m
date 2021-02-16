close all; clc; clear all;
% ElastoHydrodynamic Lubrication (EHL) Solver
% 
% The solver is based on a Finite Volume (FV) discretization of the Reynolds equation 
% incorporating mass-conserving cavitation through the cavity fraction as
% described by:
% Giacopini, M., Fowell, M. T., Dini, D., & Strozzi, A. (2010). 
% A mass-conserving complementarity formulation to study lubricant films in the presence of cavitation. 
% Journal of tribology, 132(4).
% 
% It calculates the hydrodynamic pressure and cavity fraction distribution for
% a given geometry with the Fischer-Burmeister-Newton-Schur (FBSN) algorithm
% according to:
% Woloszynski, T., Podsiadlo, P., & Stachowiak, G. W. (2015). 
% Efficient solution to the cavitation problem in hydrodynamic lubrication.
% Tribology Letters, 58(1), 18.
% The this code was largely inspired by the FBNS implementation of:
% Codrignani, A. R. (2019). Numerical representation of a pin-on-disc tribometer for the investigation of textured surfaces (Vol. 6). KIT Scientific Publishing.
%
% The code takes elastic surface deformation due to the hydrodynamic pressure into
% account through the boundary element method (BEM) applied to an elastic
% half-space model as done by:
% Pohrt, R., & Li, Q. (2014).
% Complete boundary element formulation for normal and tangential contact problems.
% Physical Mesomechanics, 17(4), 334-340.
%
% Even though this solver does not use multigrid methods, the incorporation
% of the elastic deformation into the FBNS algorithm by reducing the Kernel
% influence within the Jacobian of the FBNS algorithm to the five main
% diagonals and using the Roelands and Dowson-Higginson relations for the 
% dependence of dynamic viscosity and density of the liquid was largely inspired by:
% Venner, C. H., & Lubrecht, A. A. (Eds.). (2000). Multi-level methods in lubrication. Elsevier.
%
% Shear thinning is taken into account by limiting the Couette shear stress
% component to a maximum value by adjusting the effective dynamic viscosity
% accordingly.
%
% The load balance equation is solved using the secant algorithm. Its
% implementation in this code is a further development of the work of:
% Codrignani, A. R. (2019). Numerical representation of a pin-on-disc tribometer for the investigation of textured surfaces (Vol. 6). KIT Scientific Publishing.
%
% The initial guess of the pressure distribution is the solution of the
% dry-contact elastic half-space BEM problem. The solver is taken from the repository: 
% https://github.com/ErikHansenGit/Contact_elastic_half-space
%
% This script uses the input data specified in the script "EHL_01_setup.m" and
% exports the solution in a format which can be used by the script
% "EHL_03_visualization.m".
% The exact file path of the input data can be
% specified in the "Set file path information" section of this script. The
% file path of the output data is already specified in the script
% "EHL_01_setup.m".
% 
% Erik Hansen, 07.09.2020

%% Set file path information:
% Input path:
input_path = sprintf('%s','./../data/EHL_02_mainprocess/Input');

%% Load input information:
load(fullfile(input_path,'/alg.mat'));
load(fullfile(input_path,'/fld.mat'));
load(fullfile(input_path,'/sld.mat'));
load(fullfile(input_path,'/geo.mat'));
load(fullfile(input_path,'/opc.mat'));
load(fullfile(input_path,'/inf.mat'));
clear input_path;

% Output path:
output_main_path = inf.output_path_mainprocess;
mkdir (output_main_path)
% Save results:
output_result_path = fullfile(output_main_path,'/Result');
mkdir (output_result_path)
rmdir(output_result_path, 's')        % remove old results if they exist
mkdir (output_result_path)
% Save used input 
output_input_path = fullfile(output_main_path,'/Used_input');
mkdir (output_input_path)

maxNumCompThreads(alg.max_num_com_threads);
if alg.flag_profile_viewer
    profile on
end

% Process input data:
% Convert geometry input geo into gap height information h
h.h_g_ma    = geo.h_prof_ma;            % [m]   variation of the macroscopic gap height due to the rigid geometry
h.x1        = geo.x1;                   % [m]   coordinates in x1-direction
h.x2        = geo.x2;                   % [m]   coordinates in x2-direction
h.Nx1       = geo.Nx1;                  % [-]   number of discretization points in x1-direction
h.Nx2       = geo.Nx2;                  % [-]   number of discretization points in x1-direction
h.dx1       = geo.dx1;                  % [m]   length of discretization cell in x1-direction
h.dx2       = geo.dx2;                  % [m]   length of discretization cell in x1-direction
h.N         = h.Nx1*h.Nx2;              % [-]   total number of discretization points
h.h_el_hd_ma   = zeros(h.Nx1,h.Nx2);    % [m]   variation of the macroscopic gap height due to the elastic geometry deformation induced by hydrodynamic pressure
h.h_el_con_ma  = zeros(h.Nx1,h.Nx2);    % [m]   variation of the macroscopic gap height due to the elastic geometry deformation induced by contact pressure
h.h_pl_con_ma  = zeros(h.Nx1,h.Nx2);    % [m]   variation of the macroscopic gap height due to the plastic geometry deformation induced by contact pressure
h.h_ma      = h.h_g_ma + h.h_el_hd_ma;  % [m]   macroscopic gap height
h.h_ma_ad   = h.h_ma;                   % [m]   adjusted macroscopic gap height

% Fluid:
fld.p_W_side      = fld.p_amb;          % [Pa]  pressure at the inlet, must be in intervall of [fld.p_cav,inf] 
fld.p_E_side      = fld.p_amb;          % [Pa]  pressure at the outlet, must be in intervall of [fld.p_cav,inf] 
fld.p_S_side      = fld.p_amb;          % [Pa]  pressure on the South side, must be in intervall of [fld.p_cav,inf] 
fld.p_N_side      = fld.p_amb;          % [Pa]  pressure on the North side, must be in intervall of [fld.p_cav,inf]  

% Load:
alg.W_aim           = opc.W;            % [N]   imposed normal load force in negative x3-direction

% Save used input:
save (fullfile(output_input_path,'/fld.mat'),'fld');
save (fullfile(output_input_path,'/sld.mat'),'sld');
save (fullfile(output_input_path,'/geo.mat'),'geo');
save (fullfile(output_input_path,'/opc.mat'),'opc');
clear geo; clear output_input_path;

%% Preliminary computations and initialisations which only have to be done once:
% Start actual computation:
fprintf('\n-----------------------------\n');
fprintf('\nComputation started!\n');
fprintf('\n-----------------------------\n');

% Compute linear indices
% MATLAB first linearly indices along columns -> i -> x1 then along rows -> j -> x2:
[i_matr, j_matr]    = ndgrid(1:h.Nx1,1:h.Nx2);
h.i_lin             = zeros(h.Nx1,h.Nx2);
h.i_lin(:,:)        = sub2ind([h.Nx1 h.Nx2], i_matr(1:h.Nx1,1:h.Nx2), j_matr(1:h.Nx1,1:h.Nx2)); % [-]   linear indizes of the discretization points
clear i_matr; clear j_matr;

% Compute kernel of elastic deformations:
[h.fft2_Kernel,h.Kernel]  = construct_linear_Kernel(h.Nx1,h.Nx2,h.dx1,h.dx2,sld.E_dash);

% Initialisations:
% Stribeck curve data:
str.W               = zeros(opc.N,1);
str.F_T_low         = zeros(opc.N,1);
str.F_T_up          = zeros(opc.N,1);
str.C_f             = zeros(opc.N,1);
str.h_ma_ad_min     = zeros(opc.N,1);
str.it_tot          = zeros(opc.N,1);
str.it_W            = zeros(opc.N,1);

% Initial guess for the FBNS algorithm:
% Initialisation of the hydrodynamic pressure field:
z               = -h.h_g_ma;                        % [-]   geometry profile as an input for the contact algorithm
h_s             = alg.h_d_ma_nm2 - fld.h_min;       % [-]   seperating distance between upper and lower reference planes as an input for the contact algorithm
H               = Inf;                              % [Pa]  upper limit for the initial guess of the contact pressure
p_con_ini       = zeros(h.Nx1,h.Nx2);               % [Pa]  initial guess for contact algorithm
[p_hd_ini,~,err]= elpl_contact_pressure_akchurin_linear(fld.p_cav,H,p_con_ini,z,h_s,h.Nx1,h.Nx2,h.fft2_Kernel,alg.toll,alg.it_max_con,alg.h_ref_con);
it_p_con        = length(err);                      % [-]   number of contact algorithm iterations
res.p_con       = err;                              % [-]   relative residual of contact algorithm
clear z; clear err; clear h_s; clear H; clear p_con_ini;

%% Go through each operating condition:
for i_OC=1:opc.N   
    % Residual and algorithm data initialisation:
    res.delta_p     = zeros(alg.it_max_W*(alg.it_max_FBNS+alg.it_max_con),1);   % [-]   relative hydrodynamic pressure residual of FBNS algorithm
    res.delta_t     = zeros(alg.it_max_W*(alg.it_max_FBNS+alg.it_max_con),1);   % [-]   relative cavity fraction residual of FBNS algorithm
    res.FBNS        = zeros(alg.it_max_W*(alg.it_max_FBNS+alg.it_max_con),1);   % [-]   relative residual of FBNS algorithm
    res.W           = zeros(alg.it_max_W,2);                                    % [-]   relative residual of load balance equation
    alg.it_tot      = 1;                                                        % [-]   total iteration counter

    % Preliminary computations which only have to be done once per operating condition:
    prop.u_m        = (opc.u_low + opc.u_up(i_OC))/2;                         % [m/s] mean velocity
    prop.u_r        = opc.u_up(i_OC) - opc.u_low;                             % [m/s] relative velocity between upper and lower surface

    % Initialisation of pressures and cavity fraction:
    sol.p_hd        = p_hd_ini;                                             % [Pa]  hydrodynamic pressure
    sol.t           = zeros(h.Nx1,h.Nx2);                                   % [-]   cavity fraction
    sol.p_con       = zeros(h.Nx1,h.Nx2);                                   % [Pa]  contact pressure
    sol.p_tot       = sol.p_hd + sol.p_con;                                 % [Pa]  total pressure
    
    % Compute solutions:
    [sol,prop,h,alg,res,str] = W_bal_secant(sol,prop,h,alg,fld,str,i_OC,res);

    % Save results of each operating condition:
    sub_result_path = sprintf('/OC_%i',i_OC);
    output_sub_result_path = fullfile(output_result_path,sub_result_path);
    clear sub_result_path;
    mkdir (output_sub_result_path)
    
    save(fullfile(output_sub_result_path,'/sol.mat'),'sol');
    save(fullfile(output_sub_result_path,'/prop.mat'),'prop');
    save(fullfile(output_sub_result_path,'/h.mat'),'h');
    save(fullfile(output_sub_result_path,'/alg.mat'),'alg');
    save(fullfile(output_sub_result_path,'/res.mat'),'res');
end

% Save results of all operating conditions:
save(fullfile(output_result_path,'/str.mat'),'str');

if alg.flag_profile_viewer
    profile viewer
end

fprintf('\n-----------------------------\n');
fprintf('\nComputation finished!\n');
fprintf('\n-----------------------------\n');

%% Functions:
function [sol,prop,h,alg,res,str] = W_bal_secant(sol,prop,h,alg,fld,str,i_OC,res)
% Solve load balance equation with secant algorithm:
% The secant algorithm only suggests new rigid body displacements on the
% macroscopic scale h.h_d_ma if res.delta_p of the previous
% iteration is below alg.toll_W_delta_p and res.W is larger than alg.toll 
flag_W_bal_ini  = 0;                        % [-]   load balance initialisation flag
alg.it_W        = 1;                        % [-]   secant algorithm iteration counter
while alg.it_W == 1 || ((res.W(alg.it_W-1,2) > alg.toll || res.FBNS(alg.it_tot-1,1) > alg.toll) && alg.it_W <= alg.it_max_W) 
    % Update macroscopic rigid body displacement h.h_d_ma with secant algorithm:
    if flag_W_bal_ini == 0 % first iteration
        % Set rigid body displacement on macroscoopic scale:
        h.h_d_ma        = alg.h_d_ma_nm2;   % [m]   rigid body displacement on macroscopic scale [m]
        flag_W_bal_ini  = 1;
    elseif res.W(alg.it_W-1,2) > alg.toll && flag_W_bal_ini == 1 && res.delta_p(alg.it_tot-1,1) < alg.toll_W_delta_p % second iteration if res.r is not too large
        % Save results from previous iteration:
        alg.W_nm1       = alg.W_nm;
        alg.W_diff_nm1  = alg.W_diff_nm;
        % Set rigid body displacement on macroscopic scale:
        h.h_d_ma        = alg.h_d_ma_nm1;
        alg.h_d_ma_nm1  = alg.h_d_ma_nm2;
        flag_W_bal_ini  = 2;
    else % remaining iterations
        if res.W(alg.it_W-1,2) > alg.toll && res.delta_p(alg.it_tot-1,1) < alg.toll_W_delta_p % otherwise res.delta_p is too large
            % Save results from previous iteration:
            alg.W_nm2       = alg.W_nm1;
            alg.W_diff_nm2  = alg.W_diff_nm1;
            alg.h_d_ma_nm2  = alg.h_d_ma_nm1;
            alg.W_nm1       = alg.W_nm;
            alg.W_diff_nm1  = alg.W_diff_nm;
            alg.h_d_ma_nm1  = h.h_d_ma;
            % Compute rigid body displacement on macroscopic scale:
            h.h_d_ma        = (alg.h_d_ma_nm2*alg.W_diff_nm1 - alg.h_d_ma_nm1*alg.W_diff_nm2)/(alg.W_diff_nm1 - alg.W_diff_nm2);
        end
    end
    
    % Reset when limits are reached:
    if h.h_d_ma < alg.h_d_ma_min        % lower limit reached
        flag_W_bal_ini  = 1;
        alg.h_d_ma_nm2  = alg.h_d_ma_min;
        h.h_d_ma        = alg.h_d_ma_nm2;
        alg.h_d_ma_nm1  = alg.h_d_ma_min + alg.h_delta;
    elseif h.h_d_ma > alg.h_d_ma_max    % upper limit reached
        flag_W_bal_ini  = 1;
        alg.h_d_ma_nm2  = alg.h_d_ma_max;
        h.h_d_ma        = alg.h_d_ma_nm2;
        alg.h_d_ma_nm1  = alg.h_d_ma_max - alg.h_delta;
    end
    
    % Obtain hydrodynamic pressure field:        
    [sol,prop,h,alg,res] = FBNS(sol,prop,h,alg,res,fld);

    % Determine total pressure field:
    sol.p_tot   = sol.p_hd + sol.p_con;

    % Compute load and residual:
    alg.W_nm            = sum(sum((sol.p_tot - fld.p_amb)*h.dx1*h.dx2));        % [N]   load force in negative x3-direction
    alg.W_diff_nm       = alg.W_nm - alg.W_aim;
    res.W(alg.it_W,1)   = alg.it_tot - 1;
    res.W(alg.it_W,2)   = abs(alg.W_diff_nm/alg.W_aim);
    
    % Print iteration information:
    fprintf('\n-----------------------------\n');
    fprintf('\ni_OC      = %d\n',i_OC);
    
    fprintf('\nit_W      = %d',alg.it_W);
    fprintf('\nh_d_ma    = %d',h.h_d_ma);
    fprintf('\nW_nm      = %d',alg.W_nm);
    fprintf('\nres_W     = %d\n',res.W(alg.it_W,2));
    
    fprintf('\nit_tot    = %d',alg.it_tot-1);
    fprintf('\nres_FBNS  = %d',res.FBNS(alg.it_tot-1,1));
    fprintf('\nres_p_hd  = %d',res.delta_p(alg.it_tot-1,1));
    fprintf('\nres_t     = %d\n',res.delta_t(alg.it_tot-1,1));
    
    % Update iteration counter:
    alg.it_W = alg.it_W + 1;
end
% Save and resize final information:
alg.it_tot      = alg.it_tot - 1;
alg.it_W        = alg.it_W - 1;
res.W           = res.W(1:alg.it_W,:);
res.delta_p     = res.delta_p(1:alg.it_tot,1);
res.delta_t     = res.delta_t(1:alg.it_tot,1);
res.FBNS        = res.FBNS(1:alg.it_tot,1);

% Compute shear stresses on upper and lower surface:
[sol.tau_hd_low,sol.tau_hd_up,~] = compute_tau_hd(h,sol,prop); 
clear dp_hd_dx1;

% Stribeck data:
str.W(i_OC)             = alg.W_nm;                                 % [N]   load force in negative x3-direction
str.F_T_low(i_OC)       = sum(sum(sol.tau_hd_low*h.dx1*h.dx2));     % [N]   friction force in positive x1-direction acting from the fluid upon the lower surface 
str.F_T_up(i_OC)        = sum(sum(sol.tau_hd_up*h.dx1*h.dx2));      % [N]   friction force in positive x1-direction acting from the fluid upon the upper surface 
str.C_f(i_OC)           = str.F_T_low(i_OC)/str.W(i_OC);            % [-]   friction coefficient
str.h_ma_ad_min(i_OC)   = min(min(h.h_ma_ad));                      % [m]   minimum adjusted macroscopic gap height
str.it_tot(i_OC)        = alg.it_tot;                               % [-]   number of total iterations
str.it_W(i_OC)          = alg.it_W;                                 % [-]   number of secant iterations
end


function [sol,prop,h,alg,res] = FBNS(sol,prop,h,alg,res,fld)
% Iterative FBNS algorithm:
% For solution of the Reynolds equation with mass-conserving cavitation G(p_red,t) = A_p*p_red + B*t_red + c_G = 0 and
% Fischer-Burmeister equation F 
% p_red = p_hd - fld.p_cav reduced to a vector of length h.N
% t_red = t reduced to a vector of length h.N
% Prepare for iterations:
i           = 1:h.N;                    % [-]   counter for rows
j           = 1:h.N;                    % [-]   counter for columns
it_FBNS     = 1;                        % [-]   FBNS iteration counter
epsilon     = eps;                      % [-]   machine epsilon for F
p_red       = zeros(h.N,1); 
p_red(:)    = sol.p_hd(:) - fld.p_cav;  % [Pa] hydrodynamic pressure minus cavitation pressure in a vector in order of the linear indizes
t_red       = zeros(h.N,1);
t_red(:)    = sol.t(:);                 % [-]  cavity fraction in a vector in order of the linear indizes
while it_FBNS == 1 || (res.FBNS(alg.it_tot-1,1) > alg.toll && it_FBNS <= alg.it_max_FBNS)   
    % Compute necessary coefficients:           
    % Compute elastic gap deformation on macroscopic scale:    
    [h.h_el_hd_ma] = compute_h_el(sol.p_hd,h.Nx1,h.Nx2,h.fft2_Kernel);     
    % Update gap height on macroscopic scale
    h.h_ma      = h.h_d_ma + h.h_g_ma + h.h_el_hd_ma;     
    % Adjust gap height on macroscopic scale to prevent it from becoming zero:
    % That way, J is better conditioned and the gap height cannot become negative
    h.h_ma_ad   = h.h_ma;
    h.h_ma_ad(h.h_ma_ad<fld.h_min) = fld.h_min;     
    % Compute properties:
    prop.rho_l  = compute_rho_l(sol.p_hd,fld);                          % [kg/m^3]  density of liquid lubricant
    prop        = compute_mu_l(sol,fld,prop,h);                         % [Pas]     dynamic viscosity of liquid lubricant  
    % Compute coefficients:
    xi_p        = prop.rho_l.*h.h_ma_ad.^3./(12*prop.mu_l);             % Poiseuille coefficients for p
    xi_t        = prop.rho_l.*h.h_ma_ad*prop.u_m;                       % Couette coefficients for theta
    xi_h        = -prop.rho_l*prop.u_m;                                 % Couette coefficients for p
    
    % Construct of J_Gp, A_p, B and c_G
    [J_Gp,A_p,B,c_G] = construct_FBNS_G_matr(h,xi_p,xi_t,xi_h,fld);
    clear xi_p; clear xi_t; clear xi_h;
    
    % Compute residuals and center diagonals of J_Fp and J_Ft:
    [G,F,J_Fp_C,J_Ft_C] = compute_res(A_p,B,c_G,p_red,t_red,epsilon,h);
    clear A_p; clear c_G;
    
    % Construct the rest of J:    
    % J_Gp is reduced to the 5 main diagonals of A and J_Gt is equal to the two diagonals of B
    % Extract diagonals of J_Gp and B:
    J_Gp_S = spdiags(J_Gp,-h.Nx1);          % South diagonal
    J_Gp_W = spdiags(J_Gp,-1);              % West diagonal
    J_Gp_C = spdiags(J_Gp,0);               % Center diagonal
    J_Gp_E = spdiags(J_Gp,1);               % East diagonal
    J_Gp_N = spdiags(J_Gp,h.Nx1);           % North diagonal
    J_Gt_W = spdiags(B,-1);                 % West diagonal
    J_Gt_C = spdiags(B,0);                  % Center diagonal
    clear J_Gp; clear B;
    
    % Construct A_F, B_F, A_G and B_G:
    % Allocating arrays:
    A_F     = J_Fp_C;
    B_F     = J_Ft_C;
    A_G_S   = J_Gp_S;                     % South diagonal
    A_G_W   = J_Gp_W;                     % West diagonal
    A_G_C   = J_Gp_C;                     % Center diagonal
    A_G_E   = J_Gp_E;                     % East diagonal
    A_G_N   = J_Gp_N;                     % North diagonal
    B_G_S   = zeros(h.N,1);               % South diagonal
    B_G_W   = J_Gt_W;                     % West diagonal
    B_G_C   = J_Gt_C;                     % Center diagonal
    B_G_E   = zeros(h.N,1);               % East diagonal
    B_G_N   = zeros(h.N,1);               % North diagonal
    % Swap condition to ensure that A_F is invertible:
    sw_cond = J_Fp_C<J_Ft_C;               
    % Swap center diagonals:
    A_F(sw_cond)    = J_Ft_C(sw_cond);
    B_F(sw_cond)    = J_Fp_C(sw_cond);
    A_G_C(sw_cond)  = J_Gt_C(sw_cond);
    B_G_C(sw_cond)  = J_Gp_C(sw_cond);
    clear J_Fp_C; clear J_Ft_C; clear J_Gp_C; clear J_Gt_C;
    % Swap other diagonals:
    A_G_S(sw_cond(1:h.N-h.Nx1)) = 0;
    B_G_S(sw_cond(1:h.N-h.Nx1)) = J_Gp_S(sw_cond(1:h.N-h.Nx1));
    A_G_W(sw_cond(1:h.N-1))     = J_Gt_W(sw_cond(1:h.N-1));
    B_G_W(sw_cond(1:h.N-1))     = J_Gp_W(sw_cond(1:h.N-1));
    A_G_E(sw_cond(2:h.N))       = 0;
    B_G_E(sw_cond(2:h.N))       = J_Gp_E(sw_cond(2:h.N));
    A_G_N(sw_cond(1+h.Nx1:h.N)) = 0;
    B_G_N(sw_cond(1+h.Nx1:h.N)) = J_Gp_N(sw_cond(1+h.Nx1:h.N));
    clear J_Gp_S; clear J_Gt_W; clear J_Gp_W; clear J_Gp_E; clear J_Gp_N;
    % Construct sparse matrizes:
    A_F = sparse(i,j,A_F,h.N,h.N);
    B_F = sparse(i,j,B_F,h.N,h.N);    
    A_G = sparse(i(1+h.Nx1:h.N),j(1:h.N-h.Nx1),A_G_S(1:h.N-h.Nx1),h.N,h.N) + ...
        sparse(i(2:h.N),j(1:h.N-1),A_G_W(1:h.N-1),h.N,h.N) + ...
        sparse(i,j,A_G_C,h.N,h.N) + ...
        sparse(i(1:h.N-1),j(2:h.N),A_G_E(2:h.N),h.N,h.N) + ...
        sparse(i(1:h.N-h.Nx1),j(1+h.Nx1:h.N),A_G_N(1+h.Nx1:h.N),h.N,h.N);  
    clear A_G_S; clear A_G_W; clear A_G_C; clear A_G_E; clear A_G_N;
    B_G = sparse(i(1+h.Nx1:h.N),j(1:h.N-h.Nx1),B_G_S(1:h.N-h.Nx1),h.N,h.N) + ...
        sparse(i(2:h.N),j(1:h.N-1),B_G_W(1:h.N-1),h.N,h.N) + ...
        sparse(i,j,B_G_C,h.N,h.N) + ...
        sparse(i(1:h.N-1),j(2:h.N),B_G_E(2:h.N),h.N,h.N) + ...
        sparse(i(1:h.N-h.Nx1),j(1+h.Nx1:h.N),B_G_N(1+h.Nx1:h.N),h.N,h.N);
    clear B_G_S; clear B_G_W; clear B_G_C; clear B_G_E; clear B_G_N;
       
    % Determine delta_a and delta_b:  
    delta_b     = (B_G - A_G*(A_F\B_F))\(-G + A_G*(A_F\F));
    clear B_G; clear A_G; clear G;
    delta_a     = A_F\(-F - B_F*delta_b);
    clear A_F; clear B_F; clear F;
    
    % Determine updates for p_red and t:
    delta_p     = delta_a;
    delta_t     = delta_b;    
    % Undo swapping:
    delta_p(sw_cond) = delta_b(sw_cond);
    delta_t(sw_cond) = delta_a(sw_cond);
    clear sw_cond; clear delta_a; clear delta_b;
    
    % Update p_red and t_red:
    if it_FBNS == 1 || res.FBNS(alg.it_tot-1,1) > alg.res_alpha_FBNS
        alpha = alg.alpha_FBNS_above;
    else
        alpha = alg.alpha_FBNS_below;
    end
    p_red = p_red + alpha*delta_p;
    t_red = t_red + alpha*delta_t;
    clear alpha;
    
    % Enforce constraints:
    p_red(p_red<0) = 0;
    t_red(t_red<0) = 0;
    t_red(t_red>1) = 1;
    
    % Compute real hydrodynamic pressure field:
    sol.p_hd(:) = p_red(:) + fld.p_cav; 
    % Determine cavity fraction field:
    sol.t(:)    = t_red(:);
    
    % Compute residuals:
    res.delta_p(alg.it_tot,1) = mean(abs(delta_p))/alg.p_ref_FBNS;
    res.delta_t(alg.it_tot,1) = mean(abs(delta_t));
    clear delta_p; clear delta_t;
    res.FBNS(alg.it_tot,1)    = max([res.delta_p(alg.it_tot,1) res.delta_t(alg.it_tot,1)]);
    % Update iteration counters:
    it_FBNS     = it_FBNS + 1;
    alg.it_tot  = alg.it_tot + 1;
end
clear p_red;
clear t_red;
end

function [rho_l] = compute_rho_l(p_hd,fld)
% Dowson-Higginson relation:
% rho_l         [kg/m^3]    density of liquid lubricant
% fld.rho_0     [kg/m^3]    liquid lubricant density at ambient pressure
% p_hd          [Pa]        hydrodynamic pressure field
rho_l = fld.rho_0.*(5.9e8 + 1.34.*p_hd)./(5.9e8 + p_hd);
end

function [prop] = compute_mu_l(sol,fld,prop,h)
% Roelands relation:
% mu_l          [Pas]       dynamic viscosity of liquid lubricant
% fld.mu_0      [Pas]       liquid lubricant dynamic viscosity at ambient pressure
% p_hd          [Pa]        hydrodynamic pressure field
% fld.alpha     [1/Pa] pressure viscosity coefficient, for mineral oils between 1e-8 and 2e-8 [1/Pa]
z       = fld.alpha*1.96e8/(log(fld.mu_0) + 9.67); % [-] pressure viscosity index
prop.mu_l    = fld.mu_0.*exp((log(fld.mu_0) + 9.67)*(-1 + (1 + sol.p_hd/1.96e8).^z));
clear z;

% Consider limiting shear stress:
% Compute shear stresses at lower and upper surface and pressure gradient:
[tau_hd_low,tau_hd_up,dp_hd_dx1] = compute_tau_hd(h,sol,prop);
% Use Couette shear stress as average limit along gap height for shear shinning:
tau_av = prop.u_r./h.h_ma_ad.*prop.mu_l;
% Adjust effective dynamic viscosity:
prop.mu_l(tau_av>fld.tau_max) = fld.tau_max.*h.h_ma_ad(tau_av>fld.tau_max)/prop.u_r;
end

function [tau_hd_low,tau_hd_up,dp_hd_dx1] = compute_tau_hd(h,sol,prop)
% This function computes the shear stresses on the upper and lower surface
% and the pressure gradient:
dp_hd_dx1               = zeros(h.Nx1,h.Nx2);
i                       = 2:h.Nx1-1;
j                       = 1:h.Nx2;
dp_hd_dx1(i,j)          = (sol.p_hd(i+1,j) - sol.p_hd(i-1,j))/(2*h.dx1);                    % [Pa/m]    hydrodynamic pressure gradient in x1-direction
clear i;
dp_hd_dx1(1,j)          = (sol.p_hd(2,j) - sol.p_hd(1,j))/h.dx1;
dp_hd_dx1(h.Nx1,j)      = (sol.p_hd(h.Nx1,j) - sol.p_hd(h.Nx1-1,j))/h.dx1;
clear j;
% Assuming mu does not depend on the cavity fraction: mu = mu_l:
tau_hd_up               = (-prop.u_r./h.h_ma_ad.*prop.mu_l - 1/2*h.h_ma_ad.*dp_hd_dx1);     % [Pa]  hydrodynamic shear stress in positive x1-direction acting from the fluid upon the upper surface 
tau_hd_low              = ( prop.u_r./h.h_ma_ad.*prop.mu_l - 1/2*h.h_ma_ad.*dp_hd_dx1);     % [Pa]  hydrodynamic shear stress in positive x1-direction acting from the fluid upon the lower surface 
end

function [G,F,J_Fp_C,J_Ft_C] = compute_res(A_p,B,c_G,p_red,t_red,epsilon,h)
% This function computes the residual of the Reynolds equation G = A_p*p_red + B*t_red + c_G,
% the residual of the Fischer-Burmeister equation F = p_red + tred - sqrt(p_red^2 + t_red^2)
% and the center lines of the Jacobians J_Fp_C and J_Ft_C.
% The boundary conditions of t_red are incorporated in F, J_Fp_C and J_Ft_C and the boundary
% conditions of p_red are already incorporated in G, A_p, B and c_G
% Residual of Reynolds equation:
G = A_p*p_red + B*t_red + c_G;
% Residual of cavitation condition:
F = zeros(h.N,1);
% Inside domain excluding boundaries:
i               = 2:h.Nx1-1;           % regular indices in x1 direction
j               = 2:h.Nx2-1;           % regular indices in x2 direction
i_lin_inside    = h.i_lin(i,j);        % linear indices
clear i; clear j;
F(i_lin_inside(:)) = (p_red(i_lin_inside(:))) + (t_red(i_lin_inside(:))) ...
    - sqrt((p_red(i_lin_inside(:))).^2 + (t_red(i_lin_inside(:))).^2);
% Boundaries:
% Boundary condition for p_red is already considered in G through A_p, B and c_G
% Dirichlet of theta=0 at West boundary for theta and Neumann elsewhere is considered in F
% Boundary conditions are put up in a way that entries of J_Fp and J_Ft are >=0 or else the swapping process will not work 
% South boundary:
F(h.i_lin(1:h.Nx1,1))           = t_red(h.i_lin(1:h.Nx1,1)) - t_red(h.i_lin(1:h.Nx1,2));
% West boundary:
F(h.i_lin(1,2:h.Nx2-1)')        = t_red(h.i_lin(1,2:h.Nx2-1)');
% East boundary:
F(h.i_lin(h.Nx1,2:h.Nx2-1)')    = - t_red(h.i_lin(h.Nx1-1,2:h.Nx2-1)') + t_red(h.i_lin(h.Nx1,2:h.Nx2-1)');
% North boundary:
F(h.i_lin(1:h.Nx1,h.Nx2))       = - t_red(h.i_lin(1:h.Nx1,h.Nx2-1)) + t_red(h.i_lin(1:h.Nx1,h.Nx2));
% J_Fp and J_Ft only consist of their center diagonals of length h.N inside the domain and at the Dirichlet boundaries
% At the Theta Neumann boundaries, J_Ft is reduced to the center diagonal
% All entries of J_Fp and J_Ft must be >=0 or else the swapping process will not work 
% Residual of cavitation condition:
J_Fp_C = zeros(h.N,1);
J_Ft_C = zeros(h.N,1);
% Inside domain excluding boundaries:
J_Fp_C(i_lin_inside(:)) = 1 - p_red(i_lin_inside(:))./ ...
    sqrt((p_red(i_lin_inside(:))).^2 + (t_red(i_lin_inside(:)) + epsilon).^2);
J_Ft_C(i_lin_inside(:)) = 1 - t_red(i_lin_inside(:))./ ...
    sqrt((p_red(i_lin_inside(:))).^2 + (t_red(i_lin_inside(:)) + epsilon).^2);
clear i_lin_inside; 
% Boundaries:
% Boundary condition for p_red is already considered in G through A_p, B and c_G
% Dirichlet of theta=0 at West boundary for theta and Neumann elsewhere is considered in F 
% South boundary:
J_Ft_C(h.i_lin(1:h.Nx1,1))          = 1;
% West boundary:
J_Ft_C(h.i_lin(1,2:h.Nx2-1)')       = 1;
% East boundary:
J_Ft_C(h.i_lin(h.Nx1,2:h.Nx2-1)')   = 1;
% North boundary:
J_Ft_C(h.i_lin(1:h.Nx1,h.Nx2))      = 1;
end

function [J_Gp,A_p,B,c_G] = construct_FBNS_G_matr(h,xi_p,xi_t,xi_h,fld)
% This function conctructs the vectors and matrices needed to compute the
% residual of the Reynolds equation G = A_p*p_red + B*t_red + c_G and the
% matrix J_Gp needed to determine the update out of the residual
% The boundary conditions of t_red are later incorporated in F, J_Fp_C and J_Ft_C and the boundary
% conditions of p_red are incorporated in G, A_p, B and c_G
%
% Construction of sparse A_p:
% Discretization of Poisseuille terms with CENTRAL scheme
%
% Construction of sparse B:
% Discretization of Couette terms with UPWIND scheme
%
% Construction of c:
% Discretization of Couette terms with UPWIND scheme
%
% Construction of sparse J_Gp:
% J_Gp consists of the five main diagonals of A = A_p + A_h
% Discretization of Poisseuille terms with CENTRAL scheme and of Couette terms with UPWIND scheme
%
% Allocate arrays:
% Diagonals of A_p:
% Notation: [line, column, value] of non-sparse matrix
A_S     = zeros(h.N,3); % South diagonal
A_W     = zeros(h.N,3); % West diagonal
A_C     = zeros(h.N,3); % Center diagonal
A_E     = zeros(h.N,3); % East diagonal
A_N     = zeros(h.N,3); % North diagonal
% Diagonals of B:
% Notation: [line, column, value] of non-sparse matrix
B_W     = zeros(h.N,3); % West diagonal
B_C     = zeros(h.N,3); % Center diagonal
% Vector c:
c_G     = zeros(h.N,1);
% Coefficients:
xi_p_s  = zeros(h.Nx1,h.Nx2);
xi_p_w  = zeros(h.Nx1,h.Nx2);
xi_p_e  = zeros(h.Nx1,h.Nx2);
xi_p_n  = zeros(h.Nx1,h.Nx2);
A_s     = zeros(h.Nx1,h.Nx2);
A_w     = zeros(h.Nx1,h.Nx2);
A_e     = zeros(h.Nx1,h.Nx2);
A_n     = zeros(h.Nx1,h.Nx2);
xi_t_w  = zeros(h.Nx1,h.Nx2);
xi_t_e  = zeros(h.Nx1,h.Nx2); 
B_w     = zeros(h.Nx1,h.Nx2);
B_e     = zeros(h.Nx1,h.Nx2); 

% Inside domain excluding boundaries:
N_inside        = (h.Nx1 - 2)*(h.Nx2 - 2);  % number of points within domain
i               = 2:h.Nx1-1;                % regular indices in x1 direction
j               = 2:h.Nx2-1;                % regular indices in x2 direction
i_lin_inside    = h.i_lin(i,j);             % linear indices
% Coefficients:        
xi_p_s(i,j) = (xi_p(i,j) + xi_p(i,j-1))/2;  % south
xi_p_w(i,j) = (xi_p(i,j) + xi_p(i-1,j))/2;  % west
xi_p_e(i,j) = (xi_p(i+1,j) + xi_p(i,j))/2;  % east
xi_p_n(i,j) = (xi_p(i,j+1) + xi_p(i,j))/2;  % north
A_s(i,j)    = xi_p_s(i,j)*h.dx1/h.dx2;
A_w(i,j)    = xi_p_w(i,j)*h.dx2/h.dx1;
A_e(i,j)    = xi_p_e(i,j)*h.dx2/h.dx1;
A_n(i,j)    = xi_p_n(i,j)*h.dx1/h.dx2;
xi_t_w(i,j) = xi_t(i-1,j);                  % west
xi_t_e(i,j) = xi_t(i,j);                    % east
B_w(i,j)    = xi_t_w(i,j)*h.dx2;
B_e(i,j)    = xi_t_e(i,j)*h.dx2;
clear i; clear j;
clear xi_p_s; clear xi_p_w; clear xi_p_e; clear xi_p_n;
clear xi_u_w; clear xi_u_e; clear xi_t_w; clear xi_t_e;

% Fill diagonals of A with linear indices first along x1 then along x2:
% Notation: [line, column, value] of non-sparse matrix
A_S(1:N_inside,1:3)     = [i_lin_inside(:), i_lin_inside(:)-h.Nx1,      A_s(i_lin_inside(:))];       % South 
A_W(1:N_inside,1:3)     = [i_lin_inside(:), i_lin_inside(:)-1,          A_w(i_lin_inside(:))];       % West
A_C(1:N_inside,1:3)     = [i_lin_inside(:), i_lin_inside(:),            -A_s(i_lin_inside(:)) ...
    - A_w(i_lin_inside(:)) - A_e(i_lin_inside(:)) - A_n(i_lin_inside(:))];                           % Center
A_E(1:N_inside,1:3)     = [i_lin_inside(:), i_lin_inside(:)+1,          A_e(i_lin_inside(:))];       % East
A_N(1:N_inside,1:3)     = [i_lin_inside(:), i_lin_inside(:)+h.Nx1,      A_n(i_lin_inside(:))];       % North
% Fill diagonals of B with linear indices first along x1 then along x2: 
% Notation: [line, column, value] of non-sparse matrix
B_W(1:N_inside,1:3)     = [i_lin_inside(:), i_lin_inside(:)-1,          -B_w(i_lin_inside(:))];      % West
B_C(1:N_inside,1:3)     = [i_lin_inside(:), i_lin_inside(:),            B_e(i_lin_inside(:))];       % Center
% Fill c with linear indices first along x_1 then along x_2:
c_G(i_lin_inside(:))    = B_w(i_lin_inside(:)) - B_e(i_lin_inside(:));
clear A_s; clear A_w; clear A_e; clear A_n;
clear B_w; clear B_e;

% Boundaries:
% Dirichlet for p_red
% Coefficients for better scaling of A:
% South boundary:
k_l                     = N_inside + 1;
k_u                     = N_inside + h.Nx1;
A_C(k_l:k_u,1:3)        = [h.i_lin(1:h.Nx1,1), h.i_lin(1:h.Nx1,1), xi_p(1:h.Nx1,1)*h.dx1/h.dx2];
c_G(h.i_lin(1:h.Nx1,1)) = -xi_p(1:h.Nx1,1)*h.dx1/h.dx2*(fld.p_S_side - fld.p_cav);
% West boundary:
k_l                         = k_u + 1;
k_u                         = k_u + h.Nx2 - 2;
A_C(k_l:k_u,1:3)            = [h.i_lin(1,2:h.Nx2-1)', h.i_lin(1,2:h.Nx2-1)', xi_p(1,2:h.Nx2-1)'*h.dx2/h.dx1];
c_G(h.i_lin(1,2:h.Nx2-1)')  = -xi_p(1,2:h.Nx2-1)'*h.dx2/h.dx1*(fld.p_E_side - fld.p_cav);
% East boundary:
k_l                             = k_u + 1;
k_u                             = k_u + h.Nx2 - 2;
A_C(k_l:k_u,1:3)                = [h.i_lin(h.Nx1,2:h.Nx2-1)', h.i_lin(h.Nx1,2:h.Nx2-1)', xi_p(h.Nx1,2:h.Nx2-1)'*h.dx2/h.dx1];
c_G(h.i_lin(h.Nx1,2:h.Nx2-1)')  = -xi_p(h.Nx1,2:h.Nx2-1)'*h.dx2/h.dx1*(fld.p_W_side - fld.p_cav);
% North boundary:
k_l                         = k_u + 1;
k_u                         = k_u + h.Nx1;
A_C(k_l:k_u,1:3)            = [h.i_lin(1:h.Nx1,h.Nx2), h.i_lin(1:h.Nx1,h.Nx2), xi_p(1:h.Nx1,h.Nx2)*h.dx1/h.dx2];
c_G(h.i_lin(1:h.Nx1,h.Nx2)) = -xi_p(1:h.Nx1,h.Nx2)*h.dx1/h.dx2*(fld.p_N_side - fld.p_cav);
% Clear memory:
clear k_l;
% Communication with neighbouring cells is not needed, therefore arrays can be resized:
A_S = A_S(1:N_inside,1:3);
A_W = A_W(1:N_inside,1:3);
A_C = A_C(1:k_u,1:3);
A_E = A_E(1:N_inside,1:3);
A_N = A_N(1:N_inside,1:3);
B_W = B_W(1:N_inside,1:3);
B_C = B_C(1:N_inside,1:3);
% Clear memory:
clear k_u;

% Assemble A_p in sparse notation:
A_p = sparse(A_S(:,1),A_S(:,2),A_S(:,3),h.N,h.N) ...
    + sparse(A_W(:,1),A_W(:,2),A_W(:,3),h.N,h.N) ...
    + sparse(A_C(:,1),A_C(:,2),A_C(:,3),h.N,h.N) ...
    + sparse(A_E(:,1),A_E(:,2),A_E(:,3),h.N,h.N) ...
    + sparse(A_N(:,1),A_N(:,2),A_N(:,3),h.N,h.N);
clear A_S; clear A_W; clear A_C; clear A_E; clear A_N;

% Assemble B in sparse notation:
B = sparse(B_W(:,1),B_W(:,2),B_W(:,3),h.N,h.N) + sparse(B_C(:,1),B_C(:,2),B_C(:,3),h.N,h.N);
clear B_W; clear B_C;

% Construction of A_h with Kernel entries:
% Upwind scheme:
xi_h_w                  = zeros(h.N,1);
xi_h_e                  = zeros(h.N,1);
xi_h_w(i_lin_inside(:)) = xi_h(i_lin_inside(:)-1);
xi_h_e(i_lin_inside(:)) = xi_h(i_lin_inside(:));
% Extract Kernel entries that influence the five main diagonals of A_h:
K_S     = h.Kernel(1,2*h.Nx2);
K_SE    = h.Kernel(2,2*h.Nx2);
K_W     = h.Kernel(2*h.Nx1,1);
K_C     = h.Kernel(1,1);
K_E     = h.Kernel(2,1);
K_EE    = h.Kernel(3,1);
K_N     = h.Kernel(1,2);
K_NE    = h.Kernel(2,2);
% Allocate arrays:
h.A_h_w_aux = zeros(N_inside,3);
h.A_h_e_aux = zeros(N_inside,3);

% Fill diagonals of A_h with linear indices first along x1 then along x2:
% Notation: [line, column, value] of non-sparse matrix
A_h_S(1:N_inside,1:3)    = [i_lin_inside(:), i_lin_inside(:)-h.Nx1,       h.dx2*(-K_SE*xi_h_w(i_lin_inside(:)) + K_S*xi_h_e(i_lin_inside(:)))];      % South
A_h_W(1:N_inside,1:3)    = [i_lin_inside(:), i_lin_inside(:)-1,           h.dx2*(-K_C*xi_h_w(i_lin_inside(:)) + K_W*xi_h_e(i_lin_inside(:)))];       % West
A_h_C(1:N_inside,1:3)    = [i_lin_inside(:), i_lin_inside(:),             h.dx2*(-K_E*xi_h_w(i_lin_inside(:)) + K_C*xi_h_e(i_lin_inside(:)))];       % Center
A_h_E(1:N_inside,1:3)    = [i_lin_inside(:), i_lin_inside(:)+1,           h.dx2*(-K_EE*xi_h_w(i_lin_inside(:)) + K_E*xi_h_e(i_lin_inside(:)))];      % East
A_h_N(1:N_inside,1:3)    = [i_lin_inside(:), i_lin_inside(:)+h.Nx1,       h.dx2*(-K_NE*xi_h_w(i_lin_inside(:)) + K_N*xi_h_e(i_lin_inside(:)))];      % North
clear i_lin_inside;

% Construct sparse matrizes containing the five main diagonals:
A_h = sparse(A_h_S(:,1),A_h_S(:,2),A_h_S(:,3),h.N,h.N) + ...
    sparse(A_h_W(:,1),A_h_W(:,2),A_h_W(:,3),h.N,h.N) + ...
    sparse(A_h_C(:,1),A_h_C(:,2),A_h_C(:,3),h.N,h.N) + ...
    sparse(A_h_E(:,1),A_h_E(:,2),A_h_E(:,3),h.N,h.N) + ...
    sparse(A_h_N(:,1),A_h_N(:,2),A_h_N(:,3),h.N,h.N);
clear A_h_S; clear A_h_W; clear A_h_C; clear A_h_E; clear A_h_N; 
clear K_S; clear K_SE; clear K_W; clear K_C; clear K_E; clear K_EE; clear K_N; clear K_NE;

% Construct J_Gp:
J_Gp = A_p + A_h;
end

function [p_con,g,err] = elpl_contact_pressure_akchurin_linear(p_min,H,p_con_ini,z,h_s,Nx1,Nx2,fft2_Kernel,err_tol,it_max,h_ref)
% Erik Hansen, 26.08.2020
% Calculates the contact pressure occuring when a rigid smooth surface is
% loaded against an elastic half-space with ideal elastic-plastic material
% behaviour. The code is very similar to and consists largely of the 
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
% fft2_Kernel       [m/Pa]  2-D fast Fourier transform of the Kernel
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

function [h_el] = compute_h_el(p,Nx1,Nx2,fft2_Kernel)
% Computes elastic displacement due to pressure field with linear
% convolution in Fourier space
% Input:
% p                 [Pa]        pressure field
% Nx1               [-]         number of discretized points in x1-direction
% Nx2               [-]       	number of discretized points in x2-direction
% fft2_Kernel       [m/Pa]      2-D fast Fourier transform of the Kernel
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

function [fft2_Kernel,Kernel]  = construct_linear_Kernel(Nx1,Nx2,dx1,dx2,E_dash)
% Calculates the Kernel function for a linear convolution in the
% influence area of size Nx1*Nx2
% The Kernel is constructed for an imposed normal load as explained by Pohrt and Li, 2014:
% "Complete boundary element formulation for normal and tangential contact problems"
% under the assumption of isotropic material behaviour and therefore 
% sld.E = 2*(sld.nu + 1)*sld.G
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
% fft2_Kernel         [m/Pa]    2-D fast Fourier transform of the Kernel
% -------------------------------------------------------------------------
Nx1_K = 2*Nx1;      % [-]       Number of discretized Kernel points in x1-direction
Nx2_K = 2*Nx2;      % [-]       Number of discretized Kernel points in x2-direction
dx1_mod = dx1/2;
dx2_mod = dx2/2;
% Determine distances:
i           = 1:Nx1_K;
j           = 1:Nx2_K;
i_cond      = (i <= floor(Nx1_K/2) + 1);
x1          = (((floor(Nx1_K/2) + 1) - (i - (ceil(Nx1_K/2) + 1))) - 1)*dx1;
x1(i_cond)  = (i(i_cond) - 1)*dx1;
j_cond      = (j <= floor(Nx2_K/2) + 1);
x2          = (((floor(Nx2_K/2) + 1) - (j - (ceil(Nx2_K/2) + 1))) - 1)*dx2;
x2(j_cond)  = (j(j_cond) - 1)*dx2;
clear i; clear j; clear i_cond; clear j_cond;
clear Nx1; clear Nx2; clear dx1; clear dx2;
[x1, x2]    = ndgrid(x1,x2); % Convert vectors to matrizes
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