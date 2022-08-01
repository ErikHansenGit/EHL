close all; clc; clearvars; 
% Visualization of the EHL solver results
%
%
% Comparison to analytical solution of Fowell et al., 2007
%
% Erik Hansen, 04.03.2022
%
%==========================================================================
%% User input:
%==========================================================================

% Supply amount of simulations:
sim_N                   = 20;
% Supply simulation identifiers:
loaded_data(1).sim_id   = 'Study_E4/Nr_1/';
loaded_data(2).sim_id   = 'Study_E4/Nr_2/';
loaded_data(3).sim_id   = 'Study_E4/Nr_3/';
loaded_data(4).sim_id   = 'Study_E4/Nr_4/';
loaded_data(5).sim_id   = 'Study_E4/Nr_5/';
loaded_data(6).sim_id   = 'Study_E4/Nr_6/';
loaded_data(7).sim_id   = 'Study_E4/Nr_7/';
loaded_data(8).sim_id   = 'Study_E4/Nr_8/';
loaded_data(9).sim_id   = 'Study_E4/Nr_9/';
loaded_data(10).sim_id  = 'Study_E4/Nr_10/';
loaded_data(11).sim_id  = 'Study_E4/Nr_11/';
loaded_data(12).sim_id  = 'Study_E4/Nr_12/';
loaded_data(13).sim_id  = 'Study_E4/Nr_13/';
loaded_data(14).sim_id  = 'Study_E4/Nr_14/';
loaded_data(15).sim_id  = 'Study_E4/Nr_15/';
loaded_data(16).sim_id  = 'Study_E4/Nr_16/';
loaded_data(17).sim_id  = 'Study_E4/Nr_17/';
loaded_data(18).sim_id  = 'Study_E4/Nr_18/';
loaded_data(19).sim_id  = 'Study_E4/Nr_19/';
loaded_data(20).sim_id  = 'Study_E4/Nr_20/';

% Supply simulation descriptions:
loaded_data(1).descr    = 'UI - 3x41';
loaded_data(2).descr    = 'UI - 3x81';
loaded_data(3).descr    = 'UI - 3x161';
loaded_data(4).descr    = 'UI - 3x321';
loaded_data(5).descr    = 'UI - 3x641';
loaded_data(6).descr    = 'QUICK - 3x41';
loaded_data(7).descr    = 'QUICK - 3x81';
loaded_data(8).descr    = 'QUICK - 3x161';
loaded_data(9).descr    = 'QUICK - 3x321';
loaded_data(10).descr   = 'QUICK - 3x641';
loaded_data(11).descr   = 'UI - 3x41';
loaded_data(12).descr   = 'UI - 3x81';
loaded_data(13).descr   = 'UI - 3x161';
loaded_data(14).descr   = 'UI - 3x321';
loaded_data(15).descr   = 'UI - 3x641';
loaded_data(16).descr   = 'QUICK - 3x41';
loaded_data(17).descr   = 'QUICK - 3x81';
loaded_data(18).descr   = 'QUICK - 3x161';
loaded_data(19).descr   = 'QUICK - 3x321';
loaded_data(20).descr   = 'QUICK - 3x641';
% Choose detailed operating condition:
i_OC          = 1;

% Choose time point:
i_time        = 1;

% Output path:
output_main_path = './../data/EHL_03_visualisation/Study_E4/';

% Specify if plots should be saved and plot resolution:
flag_save_plots     =true;                                    % [-]   boolean whether to save the plots and table or not
if flag_save_plots
    print_res       = '-r600';                                  % [-]   resolution used when the plots are printed
end

% Create main output directory:
if flag_save_plots 
    warning('off','MATLAB:MKDIR:DirectoryExists')
    mkdir (output_main_path)
end

%% Load input and result properties of simulations:
for it_sim =1:sim_N
    % Main input path:
    %======================================================================
    loaded_data(it_sim).paths.input_main_path          = ...
        compose('./../data/%s/Output/', ...
        loaded_data(it_sim).sim_id);
    %======================================================================
    
    % Relative input paths:
    loaded_data(it_sim).paths.input_used_input_path    = ...
        fullfile(loaded_data(it_sim).paths.input_main_path,'Used_input/');
    loaded_data(it_sim).paths.input_result_path        = ...
        fullfile(loaded_data(it_sim).paths.input_main_path,'Result/');
    % Load used input data:
    load(fullfile(char(loaded_data(it_sim).paths.input_used_input_path),'fld.mat'));
    load(fullfile(char(loaded_data(it_sim).paths.input_used_input_path),'geo.mat'));
    load(fullfile(char(loaded_data(it_sim).paths.input_used_input_path),'opc.mat'));
    load(fullfile(char(loaded_data(it_sim).paths.input_used_input_path),'sld.mat'));
    % Save data:
    loaded_data(it_sim).fld      = fld;
    loaded_data(it_sim).geo      = geo;
    loaded_data(it_sim).opc      = opc;
    loaded_data(it_sim).sld      = sld;
    clear fld; clear geo; clear opc; clear sld; 
    % Operating condition wrapper:
    for it_OC = 1:loaded_data(it_sim).opc.N
        % Relative input paths:
        loaded_data(it_sim).loaded_OC(it_OC).paths.input_OC_path   = fullfile(char(loaded_data(it_sim).paths.input_result_path),char(compose('OC_%i/',it_OC)));
        
        % Time step wrapper:
        for it_time = 1:loaded_data(it_sim).opc.N_t(it_OC)
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).paths.input_time_path = fullfile(loaded_data(it_sim).loaded_OC(it_OC).paths.input_OC_path,char(compose('Time_%i/',it_time)));
            % Load result data:
            load(fullfile(loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).paths.input_time_path,'alg.mat'));
            load(fullfile(loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).paths.input_time_path,'h.mat'));        
            load(fullfile(loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).paths.input_time_path,'prop.mat'));
            load(fullfile(loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).paths.input_time_path,'ref.mat'));
            load(fullfile(loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).paths.input_time_path,'res.mat'));
            load(fullfile(loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).paths.input_time_path,'sol.mat'));
            % Save data:       
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).alg      = alg;
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).h        = h;
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).prop     = prop;
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).ref      = ref;
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).res      = res;
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).sol      = sol;        
            % Clear memory:
            clear alg; clear h; clear prop; clear ref; clear res; clear res; clear sol; 
        end
    end
end

%% Analytical solution:
analy.a         = loaded_data(1).geo.l_p_x1;
analy.b         = loaded_data(1).geo.l_x1;
analy.L_x1      = loaded_data(1).geo.L_x1;
analy.h_max     = loaded_data(1).geo.h_max;
analy.h_min     = loaded_data(1).geo.h_min;
analy.h_p       = loaded_data(1).geo.h_p;
analy.p_amb     = loaded_data(1).fld.p_amb;
analy.mu_l      = loaded_data(1).fld.mu_0;
analy.u_up_nc   = loaded_data(1).opc.u_up;
analy.u_up_wc   = loaded_data(11).opc.u_up;
analy.p_cav     = loaded_data(11).fld.p_cav;

analy.K         = (analy.h_max - analy.h_min)/analy.h_min;

analy.p_cf_nc   = 6*analy.u_up_nc*analy.mu_l*analy.L_x1/analy.K/analy.h_min;
analy.p_cf_wc   = 6*analy.u_up_wc*analy.mu_l*analy.L_x1/analy.K/analy.h_min;

analy.Nx1       = 40+1;
analy.x1        = linspace(0,analy.L_x1,analy.Nx1);
analy.h         = linspace(analy.h_max,analy.h_min,analy.Nx1);
analy.cond.p    = analy.x1 > analy.a     &    analy.x1 < analy.a + analy.b;
analy.h(analy.cond.p) = analy.h(analy.cond.p) + analy.h_p;
analy.h_2       = analy.h(analy.x1 == analy.a);
analy.h_2p      = analy.h_2 + analy.h_p;
analy.h_3       = analy.h(analy.x1 == analy.a + analy.b);
analy.h_3p      = analy.h_3 + analy.h_p;

analy.q_u1_nc   = (1/analy.h_2^2 - 1/analy.h_max^2) + (1/analy.h_3p^2 - 1/analy.h_2p^2) - (1/analy.h_3^2 - 1/analy.h_min^2);
analy.q_u1_nc   =((1/analy.h_2   - 1/analy.h_max)   + (1/analy.h_3p   - 1/analy.h_2p)   - (1/analy.h_3   - 1/analy.h_min)) / analy.q_u1_nc;

analy.A_1       = (analy.p_amb - analy.p_cav)/analy.p_cf_wc;

analy.q_u1_wc   = analy.h_max*analy.h_2/(analy.h_max + analy.h_2)*(analy.A_1*analy.h_max*analy.h_2/(analy.h_max - analy.h_2) + 1);

analy.A_3       = analy.A_1 + (1/analy.h_3 - 1/analy.h_min - 1/analy.h_3p) - analy.q_u1_wc*(1/analy.h_3^2 - 1/analy.h_min^2 - 1/analy.h_3p^2);
analy.h_bp      = (-1 - sqrt(1 + 4*analy.q_u1_wc*analy.A_3))/(2*analy.A_3);

analy.p_nc      = zeros(1,analy.Nx1);
analy.p_wc      = zeros(1,analy.Nx1);

analy.cond.entr = analy.x1 <= analy.a;
analy.p_nc(analy.cond.entr) = analy.p_amb + analy.p_cf_nc*((1./analy.h(analy.cond.entr) - 1/analy.h_max) - analy.q_u1_nc*(1./analy.h(analy.cond.entr).^2 - 1/analy.h_max^2));
analy.p_wc(analy.cond.entr) = analy.p_amb + analy.p_cf_wc*((1./analy.h(analy.cond.entr) - 1/analy.h_max) - analy.q_u1_wc*(1./analy.h(analy.cond.entr).^2 - 1/analy.h_max^2));

analy.p_nc_2    = analy.p_nc(analy.x1 == analy.a);
analy.p_nc(analy.cond.p) = analy.p_nc_2 + analy.p_cf_nc*((1./analy.h(analy.cond.p) - 1/analy.h_2p) - analy.q_u1_nc*(1./analy.h(analy.cond.p).^2 - 1/analy.h_2p^2));

analy.cond.cav = analy.x1 > analy.a & analy.h > analy.h_bp;
analy.p_wc(analy.cond.cav) = analy.p_cav;

analy.cond.no_cav = analy.x1 > analy.a & analy.x1 < analy.a + analy.b & analy.h <= analy.h_bp;
analy.p_wc(analy.cond.no_cav) = analy.p_cav + analy.p_cf_wc*((1./analy.h(analy.cond.no_cav) - 1/analy.h_bp) - analy.q_u1_wc*(1./analy.h(analy.cond.no_cav).^2 - 1/analy.h_bp^2));


analy.cond.exit = analy.x1 >= analy.a + analy.b;
analy.p_nc(analy.cond.exit) = analy.p_amb + analy.p_cf_nc*((1./analy.h(analy.cond.exit) - 1/analy.h_min) - analy.q_u1_nc*(1./analy.h(analy.cond.exit).^2 - 1/analy.h_min^2));
analy.p_wc(analy.cond.exit) = analy.p_amb + analy.p_cf_wc*((1./analy.h(analy.cond.exit) - 1/analy.h_min) - analy.q_u1_wc*(1./analy.h(analy.cond.exit).^2 - 1/analy.h_min^2));



% Read in data of Giacoponi:
% First comlun: x1-coordinate
% Second column: pressure of small pocket
% Third colmun: pressure of large pocket
% Woloszynski = readmatrix('../data/Woloszynski_2015/Woloszynski_2015_fig1b.csv','NumHeaderLines',1,'DecimalSeparator',',');

%==========================================================================
%% Overall plot settings:
%==========================================================================
% Set Latex style:
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

load('batlow.mat');
load('batlowS.mat');
colmap_h = batlow;
colmap_p = batlow;

batlowS_list = cell(1,100);
for i = 1:100
    batlowS_list{i} = [batlowS(i,1), batlowS(i,2), batlowS(i,3)];
end

colorlist = batlowS_list(3:10);

colors  = [1, 2, 3, 7, 8]; % These are the indices to be picked up from the colorlist
linestyle = {'-','-.','--',':'};

% Size:
widthlines          = 1.7;
contour_line_width  = 1.0;
sizeoffonts         = 9;
sizeoflegendfonts   = 9;
sizeofticksfonts    = 9;
figuresize          = [13 7];

% Ranges:
ranges.p_flag = false;
if ranges.p_flag
    ranges.p = [0 50e6]; % [Pa]
end
ranges.h_flag = true;
if ranges.h_flag
    ranges.h = [0 2e-6]; % [m]
end
ranges.x1_flag = true;
if ranges.x1_flag
    ranges.x1 = [0 2e-2]; % [m]
end
ranges.x2_flag = false;
if ranges.x2_flag
    ranges.x2 = [-6e-3 6e-3]; % [m]
end

%==========================================================================
%% Visualization:
%==========================================================================





%% Hydrodynamic pressure:
% No cavitation UI
fig = figure('Units','centimeters','Position',[21.1 11 figuresize]);
hold on
for it_sim=2:5
    
    plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.hydr_pressure_line(it_sim).p = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd(:,2);


    plot(plot_aux.hydr_pressure_line(it_sim).x1*1e3, ...
        plot_aux.hydr_pressure_line(it_sim).p*1e-6,...
        'LineWidth',widthlines,'Color',colorlist{colors(it_sim-1)},'linestyle',linestyle{it_sim-1})
        
end

plot(analy.x1*1e3,analy.p_nc*1e-6,'o','LineWidth',widthlines,'Color','k')

grid on
xlabel('$x_1 \mathrm{[mm]}$','fontsize',sizeoffonts);
ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);

ylim([0 1])
xlim([0 10])

ax = gca;
ax.FontSize = sizeofticksfonts;

legend({'$3 \cdot 81$','$3 \cdot 161$','$3 \cdot 321$','$3 \cdot 641$','Analytical'},'fontsize',sizeoffonts)

hold off  

if flag_save_plots        
    plotname = char(compose('hydr_pressure_line_nocav_UI_%i_%i',it_sim,i_OC));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;
        
% No cavitation QUICK
fig = figure('Units','centimeters','Position',[21.1 21 figuresize]);
hold on
for it_sim=7:10
    
    plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.hydr_pressure_line(it_sim).p = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd(:,2);


    plot(plot_aux.hydr_pressure_line(it_sim).x1*1e3, ...
        plot_aux.hydr_pressure_line(it_sim).p*1e-6,...
        'LineWidth',widthlines,'Color',colorlist{colors(it_sim-1-5)},'linestyle',linestyle{it_sim-1-5})
        
end

plot(analy.x1*1e3,analy.p_nc*1e-6,'o','LineWidth',widthlines,'Color','k')

grid on
xlabel('$x_1 \mathrm{[mm]}$','fontsize',sizeoffonts);
ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);

ylim([0 1])
xlim([0 10])

ax = gca;
ax.FontSize = sizeofticksfonts;

legend({'$3 \cdot 81$','$3 \cdot 161$','$3 \cdot 321$','$3 \cdot 641$','Analytical'},'fontsize',sizeoffonts)

hold off  

if flag_save_plots        
    plotname = char(compose('hydr_pressure_line_nocav_QUICK_%i_%i',it_sim,i_OC));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;
    
% With cavitation UI
fig = figure('Units','centimeters','Position',[2.1 11 figuresize]);
hold on
for it_sim=12:15
    
    plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.hydr_pressure_line(it_sim).p = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd(:,2);


    plot(plot_aux.hydr_pressure_line(it_sim).x1*1e3, ...
        plot_aux.hydr_pressure_line(it_sim).p*1e-6,...
        'LineWidth',widthlines,'Color',colorlist{colors(it_sim-1-10)},'linestyle',linestyle{it_sim-1-10})
        
end

plot(analy.x1*1e3,analy.p_wc*1e-6,'o','LineWidth',widthlines,'Color','k')

grid on
xlabel('$x_1 \mathrm{[mm]}$','fontsize',sizeoffonts);
ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);

ylim([0 10])
xlim([0 10])

ax = gca;
ax.FontSize = sizeofticksfonts;

legend({'$3 \cdot 81$','$3 \cdot 161$','$3 \cdot 321$','$3 \cdot 641$','Analytical'},'fontsize',sizeoffonts,'location','northwest')

hold off  

if flag_save_plots        
    plotname = char(compose('hydr_pressure_line_cav_UI_%i_%i',it_sim,i_OC));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;
        
% With cavitation QUICK
fig = figure('Units','centimeters','Position',[2.1 21 figuresize]);
hold on
for it_sim=17:20
    
    plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.hydr_pressure_line(it_sim).p = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd(:,2);


    plot(plot_aux.hydr_pressure_line(it_sim).x1*1e3, ...
        plot_aux.hydr_pressure_line(it_sim).p*1e-6,...
        'LineWidth',widthlines,'Color',colorlist{colors(it_sim-1-15)},'linestyle',linestyle{it_sim-1-15})
        
end

plot(analy.x1*1e3,analy.p_wc*1e-6,'o','LineWidth',widthlines,'Color','k')

grid on
xlabel('$x_1 \mathrm{[mm]}$','fontsize',sizeoffonts);
ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);

ylim([0 10])
xlim([0 10])

ax = gca;
ax.FontSize = sizeofticksfonts;

legend({'$3 \cdot 81$','$3 \cdot 161$','$3 \cdot 321$','$3 \cdot 641$','Analytical'},'fontsize',sizeoffonts,'location','northwest')

hold off  

if flag_save_plots        
    plotname = char(compose('hydr_pressure_line_cav_QUICK_%i_%i',it_sim,i_OC));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;  
    
    
    
  
% Unset Latex style:
set(groot,'defaulttextinterpreter','remove');  
set(groot,'defaultAxesTickLabelInterpreter','remove');  
set(groot,'defaultLegendInterpreter','remove');

