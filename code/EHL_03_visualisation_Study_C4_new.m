close all; clc; clearvars; 
% Visualization of the EHL solver results:
%
% Performance evaluation and comparison to Woloszynski et al., 2015
%
%
%
% Erik Hansen, 17.08.2021
%
%==========================================================================
%% User input:
%==========================================================================

% Supply amount of simulations:
sim_N                   = 32;
% Supply simulation identifiers:
loaded_data(1).sim_id   = 'Study_C4/Nr_1/';
loaded_data(2).sim_id   = 'Study_C4/Nr_2/';
loaded_data(3).sim_id   = 'Study_C4/Nr_3/';
loaded_data(4).sim_id   = 'Study_C4/Nr_4/';
loaded_data(5).sim_id   = 'Study_C4/Nr_5/';
loaded_data(6).sim_id   = 'Study_C4/Nr_6/';
loaded_data(7).sim_id   = 'Study_C4/Nr_7/';
loaded_data(8).sim_id   = 'Study_C4/Nr_8/';

loaded_data(9).sim_id   = 'Study_C4/Nr_9/';
loaded_data(10).sim_id  = 'Study_C4/Nr_10/';
loaded_data(11).sim_id  = 'Study_C4/Nr_11/';
loaded_data(12).sim_id  = 'Study_C4/Nr_12/';
loaded_data(13).sim_id  = 'Study_C4/Nr_13/';
loaded_data(14).sim_id  = 'Study_C4/Nr_14/';
loaded_data(15).sim_id  = 'Study_C4/Nr_15/';
loaded_data(16).sim_id  = 'Study_C4/Nr_16/';

loaded_data(17).sim_id  = 'Study_C4/Nr_17/';
loaded_data(18).sim_id  = 'Study_C4/Nr_18/';
loaded_data(19).sim_id  = 'Study_C4/Nr_19/';
loaded_data(20).sim_id  = 'Study_C4/Nr_20/';
loaded_data(21).sim_id  = 'Study_C4/Nr_21/';
loaded_data(22).sim_id  = 'Study_C4/Nr_22/';
loaded_data(23).sim_id  = 'Study_C4/Nr_23/';
loaded_data(24).sim_id  = 'Study_C4/Nr_24/';

loaded_data(25).sim_id  = 'Study_C4/Nr_25/';
loaded_data(26).sim_id  = 'Study_C4/Nr_26/';
loaded_data(27).sim_id  = 'Study_C4/Nr_27/';
loaded_data(28).sim_id  = 'Study_C4/Nr_28/';
loaded_data(29).sim_id  = 'Study_C4/Nr_29/';
loaded_data(30).sim_id  = 'Study_C4/Nr_30/';
loaded_data(31).sim_id  = 'Study_C4/Nr_31/';
loaded_data(32).sim_id  = 'Study_C4/Nr_32/';

% Supply simulation descriptions:
loaded_data(1).descr    = 'Ri, UI, 1';
loaded_data(2).descr    = 'Ri, UI, 2';
loaded_data(3).descr    = 'Ri, UI, 4';
loaded_data(4).descr    = 'Ri, UI, 8';
loaded_data(5).descr    = 'Ri, UI, 10';
loaded_data(6).descr    = 'Ri, UI, 20';
loaded_data(7).descr    = 'Ri, UI, 40';
loaded_data(8).descr    = 'Ri, UI, 80';

loaded_data(9).descr    = 'El, UI, 1';
loaded_data(10).descr   = 'El, UI, 2';
loaded_data(11).descr   = 'El, UI, 4';
loaded_data(12).descr   = 'El, UI, 8';
loaded_data(13).descr   = 'El, UI, 10';
loaded_data(14).descr   = 'El, UI, 20';
loaded_data(15).descr   = 'El, UI, 40';
loaded_data(16).descr   = 'El, UI, 80';

loaded_data(17).descr   = 'Ri, QUICK, 1';
loaded_data(18).descr   = 'Ri, QUICK, 2';
loaded_data(19).descr   = 'Ri, QUICK, 4';
loaded_data(20).descr   = 'Ri, QUICK, 8';
loaded_data(21).descr   = 'Ri, QUICK, 10';
loaded_data(22).descr   = 'Ri, QUICK, 20';
loaded_data(23).descr   = 'Ri, QUICK, 40';
loaded_data(24).descr   = 'Ri, QUICK, 80';

loaded_data(25).descr   = 'El, QUICK, 1';
loaded_data(26).descr   = 'El, QUICK, 2';
loaded_data(27).descr   = 'El, QUICK, 4';
loaded_data(28).descr   = 'El, QUICK, 8';
loaded_data(29).descr   = 'El, QUICK, 10';
loaded_data(30).descr   = 'El, QUICK, 20';
loaded_data(31).descr   = 'El, QUICK, 40';
loaded_data(32).descr   = 'El, QUICK, 80';

woloszynski.N = [900 3600 8100 14400 22500 32400 44100 57600 72900 90000 108900 129600 152100 176400 202500 230400 260100 291600 324900 360000];
woloszynski.t = [0.1 0.3 0.7 1.3 2.2 3.4 4.3 5.7 7.7 11.3 12.2 17.8 19.3 23.9 24.1 27.3 30.1 33.9 45.3 47.1];
% Woloszysnki only states the number of cells within the computational
% domain, while we include the boundary cells. To become consistent, the
% number of cells by woloszysnki is adjusted as follows:
woloszynski.N(:) = woloszynski.N(:) + 4*sqrt(woloszynski.N(:)) + 4;
% Choose detailed operating condition:
i_OC          = 1;

% Choose time point:
i_time        = 1;

% Output path:
output_main_path = './../data/EHL_03_visualisation/Study_C4/';

% Specify if plots should be saved and plot resolution:
flag_save_plots     = true;                                    % [-]   boolean whether to save the plots and table or not
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
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).paths.input_time_path = fullfile(loaded_data(it_sim).loaded_OC(it_OC).paths.input_OC_path,char(compose('Time_%i/',i_time)));
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

%==========================================================================
%% Overall plot settings:
%==========================================================================
% Set Latex style:
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

% batlow and batlowS Colormaps downloaded from:
% https://zenodo.org/record/4491293#.YRqGzluxVpg on 16th August 2021
%
% "The Scientific colour map batlow (Crameri 2018) is used in
% this study to prevent visual distortion of the data and exclusion of
% readers with colour­vision deficiencies (Crameri et al., 2020)."
%
% Software: Crameri, F. (2018), Scientific colour maps, Zenodo, doi:10.5281/zenodo.1243862
% Research: Crameri, F., G.E. Shephard, and P.J. Heron (2020), The mis-
% use of colour in science communication, Nature Communi-
% cations, 11, 5444. doi: 10.1038/s41467-020-19160-7

load('batlow.mat');
load('batlowS.mat');

% Surface plots:
colmap_h = batlow;
colmap_p = batlow;

% Line plots:
batlowS_list = cell(1,100);
for i = 1:100
    batlowS_list{i} = [batlowS(i,1), batlowS(i,2), batlowS(i,3)];
end
colorlist = batlowS_list(3:10);

% Size:
widthlines          = 1.7;
contour_line_width  = 1.0;
sizeoffonts         = 9;
sizeoflegendfonts   = 9;
sizeofticksfonts    = 9;
figuresize          = [13 7];

% Ranges:
ranges.p_flag = true;
if ranges.p_flag
    ranges.p = [0 0.5e6]; % [Pa]
    ticks.p = linspace(ranges.p(1),ranges.p(2),6); % [Pa]
end
ranges.h_flag = false;
if ranges.h_flag
    ranges.h = [0 5e-6]; % [m]
end
ranges.x1_flag = true;
if ranges.x1_flag
    ranges.x1 = [0e-3 8e-2]; % [m]
    ticks.x1 = linspace(ranges.x1(1),ranges.x1(2),5); % [m]
end
ranges.x2_flag = true;
if ranges.x2_flag
    ranges.x2 = [0e-3 8e-2]; % [m]
    ticks.x2 = linspace(ranges.x1(1),ranges.x2(2),5); % [m]
end

%==========================================================================
%%% Visualizations:
%==========================================================================

%% Hydrodynamic pressure fields:
for it_sim=1:sim_N   
    if it_sim == 4 
        fig = figure('Units','centimeters','Position',[01 15 figuresize]);
        subplot('Position',[0.1 0.1 0.37 0.8])
    elseif it_sim == 12 
        subplot('Position',[0.58 0.1 0.37 0.8])
    elseif it_sim == 20
        clear fig; 
        fig = figure('Units','centimeters','Position',[01 04 figuresize]);
        subplot('Position',[0.1 0.1 0.37 0.8])
    elseif it_sim == 28 
        subplot('Position',[0.58 0.1 0.37 0.8])
    end
    if it_sim == 4 || it_sim == 12 || it_sim == 20 || it_sim == 28
    
    plot_aux.hydr_pressure_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.hydr_pressure_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x2;
    plot_aux.hydr_pressure_surface(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd;
    plot_aux.hydr_pressure_surface(it_sim).p_cav = loaded_data(it_sim).fld.p_cav;
    
    surf(plot_aux.hydr_pressure_surface(it_sim).x1*1e3,plot_aux.hydr_pressure_surface(it_sim).x2*1e3,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6)
    
    hold on
    contour3(plot_aux.hydr_pressure_surface(it_sim).x1*1e3,plot_aux.hydr_pressure_surface(it_sim).x2*1e3,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6,...
        [plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6-eps plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6+eps],'LineWidth',contour_line_width,'Color',colorlist{3})        

    xlabel('$x_1 \mathrm{[mm]}$','fontsize',sizeoffonts);
    ylabel('$x_2 \mathrm{[mm]}$','fontsize',sizeoffonts);
    zlabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);
    
    shading interp
    light
    lighting 'gouraud'
    lightangle(-45,-30)
    
    colormap(colmap_p)
    ax = gca;
    ax.FontSize = sizeofticksfonts;
%     title(loaded_data(it_sim).descr + " pressure @ $U =$ " + loaded_data(it_sim).opc.u_up(i_OC) + "$\mathrm{m/s}$",'fontsize',sizeoffonts)
    if ranges.x1_flag
        xlim(ranges.x1*1e3)
        xticks(ticks.x1*1e3)
    end
    if ranges.x2_flag
        ylim(ranges.x2*1e3)
        yticks(ticks.x2*1e3)
    end
    if ranges.p_flag
        caxis(ranges.p*1e-6)
        zlim(ranges.p*1e-6)
        zticks(ticks.p*1e-6)
    end
    
    if it_sim == 4 || it_sim == 20 
        title('\textbf{(a)}')
    elseif it_sim == 12 || it_sim == 28
        title('\textbf{(b)}')
    end
    
    if flag_save_plots  
        if it_sim == 12 || it_sim == 28
            plotname = char(compose('perf_hydr_pressure_surface_%i_%i',it_sim,i_OC));
            print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
            print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
            print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
            savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
            clear plotname;
        end
    end
    clear ax;
    end
end
clear fig; 


%% EHL-FBNS Residuals:
Legend = cell(sim_N,1);
fig = figure('Units','centimeters','Position',[17 15 figuresize]);
for it_sim=1:8
    semilogy(1:loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).alg.it_tot,...
        loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.FBNS,'-','LineWidth',widthlines,'Color',colorlist{mod(it_sim-1,8)+1})
    hold on
    Legend{it_sim} = loaded_data(it_sim).descr;
end
for it_sim=9:16
    semilogy(1:loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).alg.it_tot,...
        loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.FBNS,'--','LineWidth',widthlines,'Color',colorlist{mod(it_sim-1,8)+1})
    hold on
    Legend{it_sim} = loaded_data(it_sim).descr;
end
for it_sim=17:24
    semilogy(1:loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).alg.it_tot,...
        loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.FBNS,'-.','LineWidth',widthlines,'Color',colorlist{mod(it_sim-1,8)+1})
    hold on
    Legend{it_sim} = loaded_data(it_sim).descr;
end
for it_sim=25:32
    semilogy(1:loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).alg.it_tot,...
        loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.FBNS,':','LineWidth',widthlines,'Color',colorlist{mod(it_sim-1,8)+1})
    if it_sim<sim_N
        hold on
    end
    Legend{it_sim} = loaded_data(it_sim).descr;
end
grid on
ax = gca;
ax.FontSize = sizeofticksfonts;
xlabel('$n \mathrm{[-]}$','fontsize',sizeoffonts);
ylabel('$r_{EHL-FBNS}^n \mathrm{[-]}$','fontsize',sizeoffonts);
lgd = legend(Legend);
lgd.FontSize = sizeoflegendfonts;
% if flag_save_plots 
%     plotname = 'perf_Residuals_FBNS';
%     print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
%     print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
%     print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
%     savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
%     clear plotname;
% end
clear fig; clear lgd; clear ax; clear Legend;

%% Residual definitions:
% Rigid, 1st order, 8 dimples:
% Legend = cell(9,1);
fig = figure('Units','centimeters','Position',[17 04 figuresize]);
subplot('Position',[0.125 0.3 0.35 0.6])
it_sim=4;

    semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.delta_p_nd_mean,'-','LineWidth',widthlines,'Color',colorlist{1})
    hold on

    semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.delta_p_nd_max,'-.','LineWidth',widthlines,'Color',colorlist{1})
    hold on

    semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.delta_thet_mean,'-','LineWidth',widthlines,'Color',colorlist{2})
    hold on

    semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.delta_thet_max,'-.','LineWidth',widthlines,'Color',colorlist{2})
    hold on

    semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.delta_G_max,'--','LineWidth',widthlines,'Color',colorlist{3})
    hold on

    semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.G_max,':','LineWidth',widthlines,'Color',colorlist{3})
    hold on

    semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.delta_F_max,'--','LineWidth',widthlines,'Color',colorlist{4})
    hold on
    
    semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.F_max,':','LineWidth',widthlines,'Color',colorlist{4})
    hold on

grid on
ax = gca;
ax.FontSize = sizeofticksfonts;
xlabel('$n \mathrm{[-]}$','fontsize',sizeoffonts);
ylabel('$r^n \mathrm{[-]}$','fontsize',sizeoffonts);
xlim([0 4e1])
xticks(linspace(0,4e1,5))
ylim([1e-18 1e3])
yticks([1e-18 1e-15 1e-12 1e-9 1e-6 1e-3 1e0 1e3])

title('\textbf{(a)}')


% Elastic, 1st order, 8 dimples:
Legend = cell(8,1);
subplot('Position',[0.625 0.3 0.35 0.6])
it_sim=12;
    
    dp = semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.delta_p_nd_mean,'-','LineWidth',widthlines,'Color',colorlist{1});
    hold on
    Legend{1} = '$mean,\delta p^*$';
    
    maxp = semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.delta_p_nd_max,'-.','LineWidth',widthlines,'Color',colorlist{1});
    hold on
    Legend{2} = '$max,\delta p^* $';
    
    dt = semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.delta_thet_mean,'-','LineWidth',widthlines,'Color',colorlist{2});
    hold on
    Legend{3} = '$mean,\delta \theta$';
    
    maxt = semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.delta_thet_max,'-.','LineWidth',widthlines,'Color',colorlist{2});
    hold on
    Legend{4} = '$max,\delta \theta$';
    
    dG = semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.delta_G_max,'--','LineWidth',widthlines,'Color',colorlist{3});
    hold on
    Legend{5} = '$max,\delta G$';
    
    maxG = semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.G_max,':','LineWidth',widthlines,'Color',colorlist{3});
    hold on
    Legend{6} = '$max,G$';
    
    dF = semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.delta_F_max,'--','LineWidth',widthlines,'Color',colorlist{4});
    hold on
    Legend{7} = '$max,\delta F$';
    
    maxF = semilogy(loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).res.FBNS.F_max,':','LineWidth',widthlines,'Color',colorlist{4});
    hold on
    Legend{8} = '$max,F$';

grid on
ax = gca;
ax.FontSize = sizeofticksfonts;
xlabel('$n \mathrm{[-]}$','fontsize',sizeoffonts);
ylabel('$r^n \mathrm{[-]}$','fontsize',sizeoffonts);
xlim([0 4e1])
xticks(linspace(0,4e1,5))
ylim([1e-18 1e3])
yticks([1e-18 1e-15 1e-12 1e-9 1e-6 1e-3 1e0 1e3])

title('\textbf{(b)}')

lgd = legend([dp,maxp,dt,maxt,dG,maxG,dF,maxF],...
    Legend{1},Legend{2},Legend{3},Legend{4},Legend{5},Legend{6},Legend{7},Legend{8},'Orientation','horizontal');
    set(lgd,'Position',[0.2 0.075 0.6 0.05]);
lgd.NumColumns = 4;
lgd.FontSize = sizeoflegendfonts;
    
if flag_save_plots 
    plotname =  'perf_Residuals_definition_FBNS';
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear lgd; clear ax; clear Legend;

%% Wall-clock time:
Legend = cell(5,1);
Legend{1} =  'EHL-FBNS: rigid performance';
Legend{2} =  'FBNS: original performance';
Legend{3} =  'EHL-FBNS: $N$ reference';
Legend{4} =  'EHL-FBNS: $N\log(N)$ reference';
Legend{5} =  'EHL-FBNS: $N^2$ reference';
plot_data_1 = zeros(8,2);

for it_sim=1:8
    plot_data_1(it_sim,1) = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.N;
    plot_data_1(it_sim,2) = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).alg.FBNS.eval_time.tot;
end

a = 1/( plot_data_1(1,1))^2;

b = 1/( plot_data_1(1,1));

d = 1/( plot_data_1(1,1)*log(plot_data_1(1,1)));

fig = figure('Units','centimeters','Position',[33 15 figuresize]);
loglog(plot_data_1(:,1),plot_data_1(:,2)/plot_data_1(1,2),'-','LineWidth',widthlines,'Color',colorlist{1})
hold on
loglog(woloszynski.N,woloszynski.t/woloszynski.t(1),'-','LineWidth',widthlines,'Color',colorlist{2})
hold on
loglog(plot_data_1(:,1),b*(plot_data_1(:,1)),':','LineWidth',widthlines,'Color',colorlist{3})
hold on
loglog(plot_data_1(:,1),d*(plot_data_1(:,1).*log(plot_data_1(:,1))),'-.','LineWidth',widthlines,'Color',colorlist{3})
hold on
loglog(plot_data_1(:,1),a*(plot_data_1(:,1)).^2,'--','LineWidth',widthlines,'Color',colorlist{3})

grid on
ax = gca;
ax.FontSize = sizeofticksfonts;
xlabel('$N \mathrm{[-]}$','fontsize',sizeoffonts);
ylabel('$t_{ex}^* \mathrm{[-]}$','fontsize',sizeoffonts);

xlim([1e2 1e7])
xticks([1e2 1e3 1e4 1e5 1e6 1e7])
% ylim([1e-2 1e10])
% yticks([1e-2 1e1 1e4  1e7 1e10])

lgd = legend(Legend,'Location','northwest');
lgd.FontSize = sizeoflegendfonts;
if flag_save_plots 
    plotname = 'perf_Total_Time_comp_Wol';
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear lgd; clear ax; clear Legend;


%% Wall-clock time between rigid, elastic, 1.-, 2. order and Woloszynski
% Legend = cell(5,1);
Legend = cell(4,1);
Legend{1} =  'EHL-FBNS: rigid, UI';
Legend{2} =  'EHL-FBNS: elastic, UI';
Legend{3} =  'EHL-FBNS: rigid, QUICK';
Legend{4} =  'EHL-FBNS: elastic, QUICK';
% Legend{5} =  'FBNS: Woloszynski et al.';
plot_data_1 = zeros(8,2);
plot_data_2 = zeros(8,2);
plot_data_3 = zeros(8,2);
plot_data_4 = zeros(8,2);
for it_sim=1:8
    plot_data_1(it_sim,1) = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.N;
    plot_data_1(it_sim,2) = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).alg.FBNS.eval_time.tot;
    plot_data_2(it_sim,1) = loaded_data(it_sim+8).loaded_OC(i_OC).loaded_Time(i_time).h.N;
    plot_data_2(it_sim,2) = loaded_data(it_sim+8).loaded_OC(i_OC).loaded_Time(i_time).alg.FBNS.eval_time.tot;
    plot_data_3(it_sim,1) = loaded_data(it_sim+16).loaded_OC(i_OC).loaded_Time(i_time).h.N;
    plot_data_3(it_sim,2) = loaded_data(it_sim+16).loaded_OC(i_OC).loaded_Time(i_time).alg.FBNS.eval_time.tot;
    plot_data_4(it_sim,1) = loaded_data(it_sim+24).loaded_OC(i_OC).loaded_Time(i_time).h.N;
    plot_data_4(it_sim,2) = loaded_data(it_sim+24).loaded_OC(i_OC).loaded_Time(i_time).alg.FBNS.eval_time.tot;
end


fig = figure('Units','centimeters','Position',[33 04 figuresize]);
loglog(plot_data_1(:,1),plot_data_1(:,2)/plot_data_1(1,2),'-','LineWidth',widthlines,'Color',colorlist{1})
hold on
loglog(plot_data_2(:,1),plot_data_2(:,2)/plot_data_2(1,2),':','LineWidth',widthlines,'Color',colorlist{3})
hold on
loglog(plot_data_3(:,1),plot_data_3(:,2)/plot_data_3(1,2),'--','LineWidth',widthlines,'Color',colorlist{4})
hold on
loglog(plot_data_4(:,1),plot_data_4(:,2)/plot_data_4(1,2),'-.','LineWidth',widthlines,'Color',colorlist{5})
% hold on
% loglog(woloszynski.N,woloszynski.t/woloszynski.t(1),'-','LineWidth',widthlines,'Color',colorlist{2})

grid on
ax = gca;
ax.FontSize = sizeofticksfonts;
xlabel('$N \mathrm{[-]}$','fontsize',sizeoffonts);
ylabel('$t_{ex}^* \mathrm{[-]}$','fontsize',sizeoffonts);

xlim([1e2 1e7])
xticks([1e2 1e3 1e4 1e5 1e6 1e7])
% ylim([1e-2 1e4])
% yticks([1e-2 1e-1 1e0 1e1 1e2 1e3 1e4])

lgd = legend(Legend,'Location','northwest');
lgd.FontSize = sizeoflegendfonts;
if flag_save_plots 
    plotname = 'perf_Total_Time_EHL-FBNS';
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear lgd; clear ax; clear Legend;


clear sub_result_path ;clear sub_sub_result_path; clear input_sub_sub_result_path;
% Unset Latex style:
set(groot,'defaulttextinterpreter','remove');  
set(groot,'defaultAxesTickLabelInterpreter','remove');  
set(groot,'defaultLegendInterpreter','remove');

