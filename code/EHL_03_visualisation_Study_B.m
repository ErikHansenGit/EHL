close all; clc; clearvars; 
% Visualization of the EHL solver results
%
% Rigid vs. elastic and 1st vs 2nd order evaluation and comparison to
% Bertocchi et al., 2013
%
%
%
% Erik Hansen, 04.03.2022
%
%==========================================================================
%% User input:
%==========================================================================

% Supply amount of simulations:
sim_N                   = 8;
% Supply simulation identifiers:
loaded_data(1).sim_id   = 'Study_B/Nr_1/';
loaded_data(2).sim_id   = 'Study_B/Nr_2/';
loaded_data(3).sim_id   = 'Study_B/Nr_3/';
loaded_data(4).sim_id   = 'Study_B/Nr_4/';
loaded_data(5).sim_id   = 'Study_B/Nr_5/';
loaded_data(6).sim_id   = 'Study_B/Nr_6/';
loaded_data(7).sim_id   = 'Study_B/Nr_7/';
loaded_data(8).sim_id   = 'Study_B/Nr_8/';
% Supply simulation descriptions:
loaded_data(1).descr    = 'Ri, UI, small';
loaded_data(2).descr    = 'El, UI, small';
loaded_data(3).descr    = 'Ri, UI, large';
loaded_data(4).descr    = 'El, UI, large';
loaded_data(5).descr    = 'Ri, QUICK, small';
loaded_data(6).descr    = 'El, QUICK, small';
loaded_data(7).descr    = 'Ri, QUICK, large';
loaded_data(8).descr    = 'El, QUICK, large';

% Choose detailed operating condition:
i_OC          = 1;

% Choose time point:
i_time        = 1;

% Output path:
output_main_path = './../data/EHL_03_visualisation/Study_B/';

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

% Read in data of Bertocchi:
% First comlun: x1-coordinate
% Second column: pressure of small pocket
% Third colmun: pressure of large pocket
Bertocchi = readmatrix('../data/Bertocchi_2013/Pressure_Bertocchi_2013.csv','NumHeaderLines',1,'DecimalSeparator',',');
Bertocchi(:,1) = 2e-2*Bertocchi(:,1); % convert to [m]
Bertocchi(:,2:3) = 1e6*Bertocchi(:,2:3); % convert to [Pa]

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
% ”The Scientific colour map batlow (Crameri 2018) is used in
% this study to prevent visual distortion of the data and exclusion of
% readers with colour­vision deficiencies (Crameri et al., 2020).”
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
widthlines          = 2;
contour_line_width  = 1.0;
sizeoffonts         = 9;
sizeoflegendfonts   = 9;
sizeofticksfonts    = 9;
figuresize          = [13 7];

% Ranges:
ranges.p_flag = true;
if ranges.p_flag
    ranges.p = [0 12e6]; % [Pa]
    ticks.p = linspace(ranges.p(1),ranges.p(2),7); % [Pa]
end
ranges.t_flag = true;
if ranges.t_flag
    ranges.t = [0 0.5]; % [-]
    ticks.t = linspace(ranges.t(1),ranges.t(2),6); % [-]
end
ranges.h_flag = true;
if ranges.h_flag
    ranges.h = [0 2e-6]; % [m]
    ticks.h = linspace(ranges.h(1),ranges.h(2),5); % [m]
end
ranges.x1_flag = true;
if ranges.x1_flag
    ranges.x1 = [0 2e-2]; % [m]
    ticks.x1 = linspace(ranges.x1(1),ranges.x1(2),6); % [m]
end
ranges.x2_flag = true;
if ranges.x2_flag
    ranges.x2 = [0 10e-3]; % [m]
    ticks.x2 = linspace(ranges.x1(1),ranges.x2(2),6); % [m]
end

%==========================================================================
%%% Visualization:
%==========================================================================





%% Macroscopic gap height:    
% Line plot:
for it_sim=1:sim_N
    if it_sim == 1 || it_sim == 2 || it_sim == 3 || it_sim == 5 || it_sim == 6
    
    plot_aux.h_ma_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.h_ma_surface(it_sim).Nx1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.Nx1;
    plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x2;
    plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.Nx2;
    plot_aux.h_ma_surface(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.h_ma;
    
    [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_surface(it_sim).Nx2)/2));       % [m]   x2-position of line plot
    
    end
end


%% Cavity fraction:

for it_sim=1:sim_N
    if it_sim == 1 
        fig = figure('Units','centimeters','Position',[15 15 figuresize]);
subplot('Position',[0.1 0.1 0.37 0.8])
    elseif it_sim == 2 
        subplot('Position',[0.58 0.1 0.37 0.8])
    elseif it_sim == 5
            clear fig; clear ax;
        fig = figure('Units','centimeters','Position',[15 04 figuresize]);
subplot('Position',[0.1 0.1 0.37 0.8])
    elseif it_sim == 6
        subplot('Position',[0.58 0.1 0.37 0.8])
    end
    
    
    
    if it_sim == 1 || it_sim == 2 || it_sim == 5 || it_sim == 6
    
    plot_aux.theta_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.theta_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x2;
    plot_aux.theta_surface(it_sim).theta = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.thet;
    
    surf(plot_aux.theta_surface(it_sim).x1*1e3,plot_aux.theta_surface(it_sim).x2*1e3,plot_aux.theta_surface(it_sim).theta')
    
    hold on
    plot3(plot_aux.h_ma_surface(it_sim).x1*1e3,...
        plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e3,...
        plot_aux.theta_surface(it_sim).theta(:,plot_aux.h_ma_line_x2_0_i(it_sim))','LineWidth',widthlines,'Color',colorlist{2})
    
    xlabel('$x_1 \mathrm{[mm]}$','fontsize',sizeoffonts);
    ylabel('$x_2 \mathrm{[mm]}$','fontsize',sizeoffonts);
    zlabel('$\theta \mathrm{[-]}$','fontsize',sizeoffonts);
    
    pbaspect([2 1 2])
    
    shading interp
    light
    lighting 'gouraud'
    lightangle(-45,-30)
    
    colormap(colmap_p)
    ax = gca;
    ax.FontSize = sizeofticksfonts;
    
    if it_sim == 1 || it_sim == 5 
        title('\textbf{(a)}')
    elseif it_sim == 2 || it_sim == 6
        title('\textbf{(b)}')
    end
    
    if ranges.x1_flag
        xlim(ranges.x1*1e3)
        xticks(ticks.x1*1e3)
    end
    if ranges.x2_flag
        ylim(ranges.x2*1e3)
        yticks(ticks.x2*1e3)
    end
    if ranges.t_flag
        caxis(ranges.t)
        zlim(ranges.t)
        zticks(ticks.t)
    end
    if flag_save_plots        
        plotname = char(compose('bert_theta_surface_%i_%i',it_sim,i_OC));
        print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
        print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
        print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
        savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
        clear plotname;
    end

    end

end
    clear fig; clear ax;


%% Hydrodynamic pressure:
for it_sim=1:sim_N
        if it_sim == 1 
        fig = figure('Units','centimeters','Position',[29 15 figuresize]);
subplot('Position',[0.1 0.1 0.37 0.8])
    elseif it_sim == 2 
        subplot('Position',[0.58 0.1 0.37 0.8])
    elseif it_sim == 5
            clear fig; clear ax;
        fig = figure('Units','centimeters','Position',[29 04 figuresize]);
subplot('Position',[0.1 0.1 0.37 0.8])
    elseif it_sim == 6
        subplot('Position',[0.58 0.1 0.37 0.8])
        end
    
    if it_sim == 1 || it_sim == 2 || it_sim == 5 || it_sim == 6
    
    plot_aux.hydr_pressure_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.hydr_pressure_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x2;
    plot_aux.hydr_pressure_surface(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd;
    plot_aux.hydr_pressure_surface(it_sim).p_cav = loaded_data(it_sim).fld.p_cav;
    
    surf(plot_aux.hydr_pressure_surface(it_sim).x1*1e3,plot_aux.hydr_pressure_surface(it_sim).x2*1e3,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6)
    
    hold on
    contour3(plot_aux.hydr_pressure_surface(it_sim).x1*1e3,plot_aux.hydr_pressure_surface(it_sim).x2*1e3,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6,...
        [plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6-eps plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6+eps],'LineWidth',contour_line_width,'Color',colorlist{3})        
    hold on
    plot3(plot_aux.h_ma_surface(it_sim).x1*1e3,...
        plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e3,...
        plot_aux.hydr_pressure_surface(it_sim).p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim))'*1e-6,'LineWidth',widthlines,'Color',colorlist{2})
    
    xlabel('$x_1 \mathrm{[mm]}$','fontsize',sizeoffonts);
    ylabel('$x_2 \mathrm{[mm]}$','fontsize',sizeoffonts);
    zlabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);
    
    pbaspect([2 1 2])
    
    shading interp
    light
    lighting 'gouraud'
    lightangle(-45,-30)
    
    colormap(colmap_p)
    ax = gca;
    ax.FontSize = sizeofticksfonts;
   

    if it_sim == 1 || it_sim == 5 
        title('\textbf{(a)}')
    elseif it_sim == 2 || it_sim == 6
        title('\textbf{(b)}')
    end
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
    if flag_save_plots        
        plotname = char(compose('bert_hydr_pressure_surface_%i_%i',it_sim,i_OC));
        print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
        print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
        print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
        savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
        clear plotname;
    end

    end
end
    clear fig; clear ax;



%% Comparative Line plot:
Comp_Legend = cell(4,1);
Comp_Legend{1} = 'small, EHL-FBNS';
Comp_Legend{2} = 'large, EHL-FBNS';
Comp_Legend{3} = 'small, Bertocchi et al.';
Comp_Legend{4} = 'large, Bertocchi et al.';
fig = figure('Units','centimeters','Position',[03 04 figuresize]);
it_sim = 1; 
    plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.hydr_pressure_line(it_sim).p = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim));

    
    plot(plot_aux.hydr_pressure_line(it_sim).x1*1e3, ...
        plot_aux.hydr_pressure_line(it_sim).p*1e-6,...
        'LineWidth',widthlines,'Color',colorlist{1})

    hold on
    
it_sim = 3;
    plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.hydr_pressure_line(it_sim).p = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim));

    
    plot(plot_aux.hydr_pressure_line(it_sim).x1*1e3, ...
        plot_aux.hydr_pressure_line(it_sim).p*1e-6,...
        'LineWidth',widthlines,'Color',colorlist{2})

   plot(Bertocchi(:,1)*1e3, ...
        Bertocchi(:,2)*1e-6,...
        '--','LineWidth',widthlines,'Color',colorlist{1}) 
    
    plot(Bertocchi(:,1)*1e3, ...
        Bertocchi(:,3)*1e-6,...
        '--','LineWidth',widthlines,'Color',colorlist{2}) 


    grid on
    xlabel('$x_1 \mathrm{[mm]}$','fontsize',sizeoffonts);
    ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);
    if ranges.x1_flag
        xlim(ranges.x1*1e3)
        xticks(ticks.x1*1e3)
    end
    if ranges.p_flag
        ylim([0 60])
        yticks(linspace(0,60,7))
    end
    
    ax = gca;
    ax.FontSize = sizeofticksfonts;
    
    lgd = legend(Comp_Legend{1},Comp_Legend{2},Comp_Legend{3},Comp_Legend{4},'Orientation','horizontal','Location','north');
        lgd.NumColumns = 2;
        lgd.FontSize = sizeoflegendfonts;
    
   
    if flag_save_plots        
        plotname = 'bert_hydr_pressure_line_comp_Bert';
        print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
        print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
        print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
        savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
        clear plotname;
    end
    clear fig; clear ax;
    
    
%% Line plots:

Legend = cell(4,1);
Legend{1} = 'Rigid, UI';
Legend{2} = 'Elastic, UI';
Legend{3} = 'Rigid, QUICK';
Legend{4} = 'Elastic, QUICK';

fig = figure('Units','centimeters','Position',[01 02 [figuresize(1) 3*figuresize(2)]]);

subplot('Position',[0.15 0.72 0.8 0.25])
it_sim = 1; 
    plot_aux.h_ma_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.h_ma_line(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim));

    
    ri_ui = plot(plot_aux.h_ma_line(it_sim).x1*1e3, ...
        plot_aux.h_ma_line(it_sim).h_ma*1e9,...
        'LineWidth',widthlines,'Color',colorlist{1});
    hold on
    
  it_sim = 2;
    plot_aux.h_ma_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.h_ma_line(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim));

    
    el_ui = plot(plot_aux.h_ma_line(it_sim).x1*1e3, ...
        plot_aux.h_ma_line(it_sim).h_ma*1e9,...
        'LineWidth',widthlines,'Color',colorlist{2});
    hold on
    
   it_sim = 5; 
    plot_aux.h_ma_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.h_ma_line(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim));

    
    ri_quick = plot(plot_aux.h_ma_line(it_sim).x1*1e3, ...
        plot_aux.h_ma_line(it_sim).h_ma*1e9,...
        '--','LineWidth',widthlines,'Color',colorlist{1});
    hold on
    
it_sim = 6;
    plot_aux.h_ma_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.h_ma_line(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim));

    
    el_quick = plot(plot_aux.h_ma_line(it_sim).x1*1e3, ...
        plot_aux.h_ma_line(it_sim).h_ma*1e9,...
        '--','LineWidth',widthlines,'Color',colorlist{2});
    grid on
    fig.CurrentAxes.YDir = 'Reverse';
    xlabel('$x_1 \mathrm{[mm]}$','fontsize',sizeoffonts);
    ylabel('$h \mathrm{[nm]}$','fontsize',sizeoffonts);
    if ranges.x1_flag
        xlim(ranges.x1*1e3)
        xticks(ticks.x1*1e3)
    end
    if ranges.h_flag
        ylim(ranges.h*1e9)
        yticks(ticks.h*1e9)
    end
    ax = gca;
    ax.FontSize = sizeofticksfonts;
    
    title('\textbf{(a)}')
   

subplot('Position',[0.15 0.385 0.8 0.25])
 it_sim = 1;
        plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
        plot_aux.hydr_pressure_line(it_sim).p = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim));


        plot(plot_aux.hydr_pressure_line(it_sim).x1*1e3, ...
            plot_aux.hydr_pressure_line(it_sim).p*1e-6,...
            'LineWidth',widthlines,'Color',colorlist{1})
        hold on
        
       it_sim = 2; 
        plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
        plot_aux.hydr_pressure_line(it_sim).p = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim));

        plot(plot_aux.hydr_pressure_line(it_sim).x1*1e3, ...
            plot_aux.hydr_pressure_line(it_sim).p*1e-6,...
            'LineWidth',widthlines,'Color',colorlist{2})
        hold on
        
        it_sim = 5;
        plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
        plot_aux.hydr_pressure_line(it_sim).p = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim));

        plot(plot_aux.hydr_pressure_line(it_sim).x1*1e3, ...
            plot_aux.hydr_pressure_line(it_sim).p*1e-6,...
            '--','LineWidth',widthlines,'Color',colorlist{1})
        hold on
        
        it_sim = 6;
        plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
        plot_aux.hydr_pressure_line(it_sim).p = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim));

        
        plot(plot_aux.hydr_pressure_line(it_sim).x1*1e3, ...
            plot_aux.hydr_pressure_line(it_sim).p*1e-6,...
            '--','LineWidth',widthlines,'Color',colorlist{2})
        
        
        grid on
        xlabel('$x_1 \mathrm{[mm]}$','fontsize',sizeoffonts);
        ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);
        if ranges.x1_flag
            xlim(ranges.x1*1e3)
            xticks(ticks.x1*1e3)
        end
        if ranges.p_flag
            ylim(ranges.p*1e-6)
            yticks(ticks.p*1e-6)
        end

        ax = gca;
        ax.FontSize = sizeofticksfonts;
        
        title('\textbf{(b)}')
        

subplot('Position',[0.15 0.06 0.8 0.25])
it_sim = 1;
        plot_aux.theta_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
        plot_aux.theta_line(it_sim).theta = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.thet(:,plot_aux.h_ma_line_x2_0_i(it_sim));

        plot(plot_aux.theta_line(it_sim).x1*1e3, ...
            plot_aux.theta_line(it_sim).theta,...
            'LineWidth',widthlines,'Color',colorlist{1})
            hold on
            
            it_sim = 2;
        plot_aux.theta_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
        plot_aux.theta_line(it_sim).theta = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.thet(:,plot_aux.h_ma_line_x2_0_i(it_sim));

        plot(plot_aux.theta_line(it_sim).x1*1e3, ...
            plot_aux.theta_line(it_sim).theta,...
            'LineWidth',widthlines,'Color',colorlist{2})
            hold on
            
           it_sim = 5;
        plot_aux.theta_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
        plot_aux.theta_line(it_sim).theta = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.thet(:,plot_aux.h_ma_line_x2_0_i(it_sim));

        plot(plot_aux.theta_line(it_sim).x1*1e3, ...
            plot_aux.theta_line(it_sim).theta,...
            '--','LineWidth',widthlines,'Color',colorlist{1})
            hold on
            
          it_sim = 6;
        plot_aux.theta_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
        plot_aux.theta_line(it_sim).theta = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.thet(:,plot_aux.h_ma_line_x2_0_i(it_sim));

        plot(plot_aux.theta_line(it_sim).x1*1e3, ...
            plot_aux.theta_line(it_sim).theta,...
            '--','LineWidth',widthlines,'Color',colorlist{2})


 grid on
        xlabel('$x_1 \mathrm{[mm]}$','fontsize',sizeoffonts);
        ylabel('$\theta \mathrm{[-]}$','fontsize',sizeoffonts);
        if ranges.x1_flag
            xlim(ranges.x1*1e3)
            xticks(ticks.x1*1e3)
        end

        if ranges.t_flag
            ylim(ranges.t)
            yticks(ticks.t)
        end
        ax = gca;
        ax.FontSize = sizeofticksfonts;
        
        title('\textbf{(c)}')
        
        lgd = legend(Legend{1},Legend{2},Legend{3},Legend{4},'Orientation','horizontal','Location','north');
        lgd.NumColumns = 2;
        lgd.FontSize = sizeoflegendfonts;
        set(lgd,'Position',[0.25 0.915 0.6 0.05]);
    
    
        if flag_save_plots        
            plotname = 'bert_line_B';
            print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
            print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
            print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
            savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
            clear plotname;
        end
        clear fig; clear ax;

clear sub_result_path ;clear sub_sub_result_path; clear input_sub_sub_result_path;
% Unset Latex style:
set(groot,'defaulttextinterpreter','remove');  
set(groot,'defaultAxesTickLabelInterpreter','remove');  
set(groot,'defaultLegendInterpreter','remove');

