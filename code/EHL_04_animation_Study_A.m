close all; clc; clear all; 
% Video maker of the EHL solver results
%
% Videos are exported as .avi files
%
% Altay Kacan, Erik Hansen 
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Execution of this script takes quite long because a lot of data has to be
% read in!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%==========================================================================
%% User input:
%==========================================================================

% Output path:
output_main_path ='./../data/EHL_03_visualisation/Study_A/Videos';

warning('off','MATLAB:MKDIR:DirectoryExists')
mkdir (output_main_path)
    
% Supply the amount of total simulations:
sim_N                   = 8;

% Supply simulation identifiers:
loaded_data(1).sim_id   = 'Study_A/Nr_1/';
loaded_data(2).sim_id   = 'Study_A/Nr_2/';
loaded_data(3).sim_id   = 'Study_A/Nr_3/';
loaded_data(4).sim_id   = 'Study_A/Nr_4/';
loaded_data(5).sim_id   = 'Study_A/Nr_5/';
loaded_data(6).sim_id   = 'Study_A/Nr_6/';
loaded_data(7).sim_id   = 'Study_A/Nr_7/';
loaded_data(8).sim_id   = 'Study_A/Nr_8/';

% Supply simulation descriptions:
loaded_data(1).descr    = 'Deep, UI, SSR = 0';
loaded_data(2).descr    = 'Shallow, UI, SSR = 0';
loaded_data(3).descr    = 'Deep, UI, SSR = -0.5';
loaded_data(4).descr    = 'Shallow, UI, SSR = -0.5';
loaded_data(5).descr    = 'Deep, QUICK, SSR = 0';
loaded_data(6).descr    = 'Shallow, QUICK, SSR = 0';
loaded_data(7).descr    = 'Deep, QUICK, SSR = -0.5';
loaded_data(8).descr    = 'Shallow, QUICK, SSR = -0.5';

% Choose which simulations to create a video for:
sims_to_plot = [1 2 3 4 5 6 7 8];

% Choose detailed operating condition:
i_OC          = 1;

% Define the frames per second of the video(FPS):
fps = 30;

% Define flags for video settings:
flag_vid.h_surf       = true; % gap height distribution for the whole contact
flag_vid.p_surf       = true; % pressure distribution for the whole contact
flag_vid.extra_info   = true; % minimum gap height and maximum pressure
flag_vid.load_force   = false; % load force, doesn't make sense to plot since it stays constant due to load balance equation

ranges.h    = [0 12];       % [mu m]    range for p-axis and colormap in gap height surface plots
ranges.p_ax = [0 1500];     % [MPa] range for p-axis in pressure surface plots 
ranges.p_co = [0 500];      % [MPa] range for colormap in pressure surface plots 

%==========================================================================
%% Loading the data:
%==========================================================================    
for it_sim = 1:sim_N
   if ismember(it_sim, sims_to_plot)
       % Main input path:
       %==========================================================================
        loaded_data(it_sim).paths.input_main_path          = ...
            compose('./../data/%s/Output/', ...
            loaded_data(it_sim).sim_id);
       %==========================================================================
       % Define relative input paths (for used input & result):
        loaded_data(it_sim).paths.input_used_input_path    = ...
            fullfile(loaded_data(it_sim).paths.input_main_path,'Used_input/');
        loaded_data(it_sim).paths.input_result_path        = ...
            fullfile(loaded_data(it_sim).paths.input_main_path,'Result/');

        % Load used input data:
        load(fullfile(char(loaded_data(it_sim).paths.input_used_input_path),'fld.mat'));
        load(fullfile(char(loaded_data(it_sim).paths.input_used_input_path),'geo.mat'));
        load(fullfile(char(loaded_data(it_sim).paths.input_used_input_path),'opc.mat'));
        load(fullfile(char(loaded_data(it_sim).paths.input_used_input_path),'sld.mat'));        

        % Save data the loaded data:
        loaded_data(it_sim).fld      = fld;
        loaded_data(it_sim).geo      = geo;
        loaded_data(it_sim).opc      = opc;
        loaded_data(it_sim).sld      = sld;        
        clear fld; clear geo; clear opc; clear sld; 

        % Operating condition wrapper - goes through all operating conditions:
        for it_OC = 1:loaded_data(it_sim).opc.N
           % Relative input paths for this operating condition:
           loaded_data(it_sim).loaded_OC(it_OC).paths.input_OC_path  = ...
               fullfile(char(loaded_data(it_sim).paths.input_result_path),char(compose('OC_%i/',it_OC))); 

           % Time step wrapper - goes through all time points:
           for it_time = 1:loaded_data(it_sim).opc.N_t(it_OC) 
%            for it_time  = 1:15 %TEST
           fprintf('\nRead in-------------------------');
           fprintf('\nit_sim        = %i',it_sim);
           fprintf('\nit_OC         = %i',it_OC);
           fprintf('\nit_time       = %i',it_time);
           fprintf('\n--------------------------------');
            % Define the paths: 
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).paths.input_time_path = ...
                fullfile(loaded_data(it_sim).loaded_OC(it_OC).paths.input_OC_path,char(compose('Time_%i/',it_time)));
            % Load used input data:
            load(fullfile(char(loaded_data(it_sim).paths.input_used_input_path),'h_time',sprintf('OC_%i',it_OC),sprintf('Time_%i',it_time),'slc.mat'));
            % Load result data:
            load(fullfile(loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).paths.input_time_path,'alg.mat'));
            load(fullfile(loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).paths.input_time_path,'h.mat'));        
            load(fullfile(loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).paths.input_time_path,'prop.mat'));
            load(fullfile(loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).paths.input_time_path,'ref.mat'));
            load(fullfile(loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).paths.input_time_path,'res.mat'));
            load(fullfile(loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).paths.input_time_path,'sol.mat'));

            % Save the loaded data:       
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).alg      = alg;
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).h        = h;
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).prop     = prop;
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).ref      = ref;
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).res      = res;
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).sol      = sol;                  
            loaded_data(it_sim).loaded_OC(it_OC).loaded_Time(it_time).slc      = slc; 
            % Clear memory:
            clear alg; clear h; clear prop; clear ref; clear res; clear res; clear sol; clear slc;                
           end
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

% % Line plots:
% KIT_colorlist={[0,150,130]/255,[162 34 35]/255,[70 100 170]/255,[252 229 0]/255,[140 182 60]/256,[223 155 27]/255,[167 130 46]/255,[163 16 124]/255,[35 161 224]/255};
% % Surface plots:
% colmap_h = parula;
% colmap_p = parula;

% batlow and batlowS Colormaps downloaded from:
% https://zenodo.org/record/4491293#.YRqGzluxVpg on 16th August 2021
%
% â€?The Scientific colour map batlow (Crameri 2018) is used in
% this study to prevent visual distortion of the data and exclusion of
% readers with colourÂ­vision deficiencies (Crameri et al., 2020).â€?
%
% Software: Crameri, F. (2018), Scientific colour maps, Zenodo, doi:10.5281/zenodo.1243862
% Research: Crameri, F., G.E. Shephard, and P.J. Heron (2020), The mis-
% use of colour in science communication, Nature Communi-

load('batlow.mat');
load('batlowS.mat');
colmap_h = batlow;
colmap_p = batlow;

batlowS_list = cell(1,100);
for i = 1:100
    batlowS_list{i} = [batlowS(i,1), batlowS(i,2), batlowS(i,3)];
end

KIT_colorlist = batlowS_list(3:10);

% Size:
widthlines          = 2;
contour_line_width  = 1.0;
sizeoffonts         = 18;
sizeoflegendfonts   = 18;
sizeofticksfonts    = 18;
figuresize          = [8 5];

%% Figure Creation:
%
% config 1 : all plots
% config 2 : both surface plots with extra info
% config 3 : both surface plots with load force (no extra info)
% config 4 : both surface plots only
% config 5 : load force only
% config 6 : separate surface plot only
% config 7 : separate surface plot with extra info

if flag_vid.h_surf && flag_vid.p_surf && flag_vid.extra_info && flag_vid.load_force
    config  = 1; 
    rows    = 3;
    columns = 2;
elseif flag_vid.h_surf && flag_vid.p_surf && flag_vid.extra_info && (flag_vid.load_force == false)
    config  = 2;
    rows    = 2;
    columns = 2;    
elseif flag_vid.h_surf && flag_vid.p_surf && flag_vid.load_force && (flag_vid.extra_info == false)
    config  = 3;
    rows    = 2;
    columns = 2; 
elseif flag_vid.h_surf && flag_vid.p_surf && ((flag_vid.load_force || flag_vid.extra_info) == false)
    config  = 4;
    rows    = 1;
    columns = 2;    
elseif flag_vid.load_force && (( flag_vid.extra_info || flag_vid.h_surf || flag_vid.p_surf) == false)
    config  = 5;
    rows    = 1;
    columns = 2;
elseif (flag_vid.h_surf || flag_vid.p_surf) && ((flag_vid.load_force || flag_vid.extra_info ) == false)
    config  = 6;
    rows    = 1;
    columns = 1;    
elseif (flag_vid.h_surf || flag_vid.p_surf) && flag_vid.extra_info && (flag_vid.load_force == false)
    config  = 7;
    rows    = 2;
    columns = 1;
end


fprintf("configuration: " + config + " \n")

%% Plotting and Video Creation: 
for it_sim = 1:sim_N
   if ismember(it_sim, sims_to_plot)
       % Create figure based on plots to show:
    
       fig  = figure('Units','normalized','Position', [0 0 1 1]); 
       
       % Get how many time points are needed:
       time_N = loaded_data(it_sim).opc.N_t(it_OC); 
%        time_N  = 15; %TEST
       % Create the necessary variables for the minimum gap height, maximum
       % pressure and load force graphs:
    
       for it_time = 1:time_N
           h_for_this_time     = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).h.h_ma;
           p_for_this_time     = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).sol.p_hd;
           force_for_this_time = sum(p_for_this_time*loaded_data(it_sim).geo.dx1*loaded_data(it_sim).geo.dx2 ...
               ,'all'); % [N] the generated load carrying force is the (discrete) integral of the pressure field
           goal_force_for_this_time = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).sol.F_N;  
           x1_c_for_this_time = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).slc.x1_c;  
           plot_aux.min_h(it_sim).value(it_time)         = min(min(h_for_this_time));
           plot_aux.min_h(it_sim).dimple_center(it_time) = x1_c_for_this_time;       
           plot_aux.max_p(it_sim).value(it_time)         = max(max(p_for_this_time));
           plot_aux.max_p(it_sim).dimple_center(it_time) = x1_c_for_this_time; 
           plot_aux.load(it_sim).value(it_time)          = force_for_this_time;
           plot_aux.load(it_sim).goal(it_time)           = goal_force_for_this_time;
           plot_aux.load(it_sim).dimple_center(it_time)  = x1_c_for_this_time;
       end
       clear h_for_this_time; clear p_for_this_time; clear force_for_this_time; clear goal_force_for_this_time;
       clear x1_c_for_this_time;
       % Set the time resolution for the video:
       time_reso = time_N;
       % Start time loop for plotting and animating:
       for it_time = 1:time_N
           % Clear the figure:
           clf
           % Plot the required subplots depending on the configuration:
           switch config
               case 1
                  %========================================================================== 
                  h_surf = subplot(rows, columns, 1);
                   plot_aux.h_ma_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x1;
                   plot_aux.h_ma_surface(it_sim).Nx1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.Nx1;
                   plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x2;
                   plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.Nx2;
                   plot_aux.h_ma_surface(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.h_ma;       

                   %[plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_surface(it_sim).Nx2)/2));       % [m]   x2-position of line plot
                   [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(ceil(plot_aux.h_ma_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

                   surf(plot_aux.h_ma_surface(it_sim).x1*1e6,plot_aux.h_ma_surface(it_sim).x2*1e6,plot_aux.h_ma_surface(it_sim).h_ma'*1e6)
                   hold on
                   plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
                    plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
                    plot_aux.h_ma_surface(it_sim).h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim))'*1e6,'LineWidth',widthlines,'Color',KIT_colorlist{2})

                   xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   zlabel('$h \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   xticks(-500:250:500);
                   yticks(-500:250:500);
                   zlim(ranges.h);
                   fig.CurrentAxes.ZDir = 'Reverse';
                   shading interp
                   light
                   lighting 'gouraud'
                   lightangle(-45,-30)
                   colormap(h_surf,flipud(colmap_h))
                   caxis(ranges.h)
                   ax = gca;
                   ax.FontSize = sizeofticksfonts;
                   title("$h$");
                  
                  %==========================================================================
                  p_surf = subplot(rows, columns, 2);
                   plot_aux.hydr_pressure_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x1;
                   plot_aux.hydr_pressure_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x2;
                   plot_aux.hydr_pressure_surface(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).sol.p_hd;
                   plot_aux.hydr_pressure_surface(it_sim).p_cav = loaded_data(it_sim).fld.p_cav;
                   
                   surf(plot_aux.hydr_pressure_surface(it_sim).x1*1e6,plot_aux.hydr_pressure_surface(it_sim).x2*1e6,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6)

                   hold on
                   contour3(plot_aux.hydr_pressure_surface(it_sim).x1*1e6,plot_aux.hydr_pressure_surface(it_sim).x2*1e6,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6,...
                       [plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6-eps plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6+eps],'LineWidth',widthlines,'Color',KIT_colorlist{3})        
                   hold on
                   plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
                       plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
                       plot_aux.hydr_pressure_surface(it_sim).p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim))'*1e-6,'LineWidth',widthlines,'Color',KIT_colorlist{2})

                   xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   zlabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);
                   xticks(-500:250:500);
                   yticks(-500:250:500);  
                   xlim([-500 500]);
                   ylim([-500 500]);
                   zlim(ranges.p_ax);
                   shading interp
                   light
                   lighting 'gouraud'
                   lightangle(-45,-30)

                   colormap(p_surf,colmap_p)
                   caxis(ranges.p_co)
                   ax = gca;
                   ax.FontSize = sizeofticksfonts;                  
                   title("$p_{hd}$");

                  %==========================================================================
                  h_min = subplot(rows, columns, 3);
                   plot(plot_aux.min_h(it_sim).dimple_center*1e6, plot_aux.min_h(it_sim).value*1e9, 'LineWidth',widthlines,'Color',KIT_colorlist{1} );
                   hold on
                   grid on
                   % Red dot:
                   plot(plot_aux.min_h(it_sim).dimple_center(it_time)*1e6, plot_aux.min_h(it_sim).value(it_time)*1e9,'ro','MarkerFaceColor','red','MarkerSize',2);                 
                   % Vertical line:
                   line([plot_aux.min_h(it_sim).dimple_center(it_time)*1e6 plot_aux.min_h(it_sim).dimple_center(it_time)*1e6],[0 plot_aux.min_h(it_sim).value(it_time)*1e9]...
                       ,'LineStyle',':','LineWidth',widthlines/2,'Color','red');                   
                   xlabel('$x_{1,c}$ [$\mu$m]','fontsize',sizeoffonts);
                   ylabel('$h_{min}$ [nm]','fontsize',sizeoffonts);   
                   xticks(-500:250:500);
                   xlim([-500 500]);
                   ylim([0 300]);
                   
                   ax = gca;
                   ax.FontSize = sizeofticksfonts;
                   title("$h_{min}$");

                  %==========================================================================
                  p_max = subplot(rows, columns, 4);
                   plot(plot_aux.max_p(it_sim).dimple_center*1e6, plot_aux.max_p(it_sim).value*1e-6, 'LineWidth',widthlines,'Color',KIT_colorlist{1} );
                   hold on
                   grid on
                   % Red dot:
                   plot(plot_aux.max_p(it_sim).dimple_center(it_time)*1e6, plot_aux.max_p(it_sim).value(it_time)*1e-6,'ro','MarkerFaceColor','red','MarkerSize',2);
                   % Vertical line:
                   line([plot_aux.max_p(it_sim).dimple_center(it_time)*1e6 plot_aux.max_p(it_sim).dimple_center(it_time)*1e6],[0 plot_aux.max_p(it_sim).value(it_time)*1e-6]...
                       ,'LineStyle',':','LineWidth',widthlines/2,'Color','red');
                   
                   xlabel('$x_{1,c}$ [$\mu$m]','fontsize',sizeoffonts);
                   ylabel('$p_{hd,max}$ [MPa]','fontsize',sizeoffonts);   
                   xticks(-500:250:500);
                   yticks(0:250:1500);
                   xlim([-500 500]);
                   ylim([0 1500]);

                   ax = gca;
                   ax.FontSize = sizeofticksfonts;
                   title("$p_{hd, max}$");

                  %==========================================================================
                  load_force = subplot(rows,columns, [5,6]);
                   % I wasn't sure what to plot here exactly because it is
                   % almost constant
                   title("Load Force");

                  %==========================================================================
               case 2
                  %========================================================================== 
                  h_surf = subplot(rows, columns, 1);
                   plot_aux.h_ma_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x1;
                   plot_aux.h_ma_surface(it_sim).Nx1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.Nx1;
                   plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x2;
                   plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.Nx2;
                   plot_aux.h_ma_surface(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.h_ma;       

                   %[plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_surface(it_sim).Nx2)/2));       % [m]   x2-position of line plot
                   [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(ceil(plot_aux.h_ma_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

                   surf(plot_aux.h_ma_surface(it_sim).x1*1e6,plot_aux.h_ma_surface(it_sim).x2*1e6,plot_aux.h_ma_surface(it_sim).h_ma'*1e6)
                   hold on
                   plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
                    plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
                    plot_aux.h_ma_surface(it_sim).h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim))'*1e6,'LineWidth',widthlines,'Color',KIT_colorlist{2})

                   xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   zlabel('$h \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   xticks(-500:250:500);
                   yticks(-500:250:500);
                   zlim(ranges.h);
                   fig.CurrentAxes.ZDir = 'Reverse';
                   shading interp
                   light
                   lighting 'gouraud'
                   lightangle(-45,-30)
                   colormap(h_surf,flipud(colmap_h))
                   caxis(ranges.h)
                   ax = gca;
                   ax.FontSize = sizeofticksfonts;
                   title("$h$");
                  %==========================================================================
                  p_surf = subplot(rows, columns, 2);
                   plot_aux.hydr_pressure_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x1;
                   plot_aux.hydr_pressure_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x2;
                   plot_aux.hydr_pressure_surface(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).sol.p_hd;
                   plot_aux.hydr_pressure_surface(it_sim).p_cav = loaded_data(it_sim).fld.p_cav;
                   
                   surf(plot_aux.hydr_pressure_surface(it_sim).x1*1e6,plot_aux.hydr_pressure_surface(it_sim).x2*1e6,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6)

                   hold on
                   contour3(plot_aux.hydr_pressure_surface(it_sim).x1*1e6,plot_aux.hydr_pressure_surface(it_sim).x2*1e6,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6,...
                       [plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6-eps plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6+eps],'LineWidth',widthlines,'Color',KIT_colorlist{3})        
                   hold on
                   plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
                       plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
                       plot_aux.hydr_pressure_surface(it_sim).p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim))'*1e-6,'LineWidth',widthlines,'Color',KIT_colorlist{2})

                   xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   zlabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);
                   xticks(-500:250:500);
                   yticks(-500:250:500);  
                   xlim([-500 500]);
                   ylim([-500 500]);
                   zlim(ranges.p_ax);
                   shading interp
                   light
                   lighting 'gouraud'
                   lightangle(-45,-30)

                   colormap(p_surf,colmap_p)
                   caxis(ranges.p_co)
                   ax = gca;
                   ax.FontSize = sizeofticksfonts;                  
                   title("$p_{hd}$");                   
                   
                  %==========================================================================
                  h_min = subplot(rows, columns, 3);
                   plot(plot_aux.min_h(it_sim).dimple_center*1e6, plot_aux.min_h(it_sim).value*1e9, 'LineWidth',widthlines,'Color',KIT_colorlist{1} );
                   hold on
                   grid on
                   % Red dot:
                   plot(plot_aux.min_h(it_sim).dimple_center(it_time)*1e6, plot_aux.min_h(it_sim).value(it_time)*1e9,'ro','MarkerFaceColor','red','MarkerSize',2);                 
                   % Vertical line:
                   line([plot_aux.min_h(it_sim).dimple_center(it_time)*1e6 plot_aux.min_h(it_sim).dimple_center(it_time)*1e6],[0 plot_aux.min_h(it_sim).value(it_time)*1e9]...
                       ,'LineStyle',':','LineWidth',widthlines/2,'Color','red');                   
                   xlabel('$x_{1,c}$ [$\mu$m]','fontsize',sizeoffonts);
                   ylabel('$h_{min}$ [nm]','fontsize',sizeoffonts);   
                   xticks(-500:250:500);
                   xlim([-500 500]);
                   ylim([0 300]);
                   
                   ax = gca;
                   ax.FontSize = sizeofticksfonts;
                   title("$h_{min}$");
                  %==========================================================================
                  p_max = subplot(rows, columns, 4);
                   plot(plot_aux.max_p(it_sim).dimple_center*1e6, plot_aux.max_p(it_sim).value*1e-6, 'LineWidth',widthlines,'Color',KIT_colorlist{1} );
                   hold on
                   grid on
                   % Red dot:
                   plot(plot_aux.max_p(it_sim).dimple_center(it_time)*1e6, plot_aux.max_p(it_sim).value(it_time)*1e-6,'ro','MarkerFaceColor','red','MarkerSize',2);
                   % Vertical line:
                   line([plot_aux.max_p(it_sim).dimple_center(it_time)*1e6 plot_aux.max_p(it_sim).dimple_center(it_time)*1e6],[0 plot_aux.max_p(it_sim).value(it_time)*1e-6]...
                       ,'LineStyle',':','LineWidth',widthlines/2,'Color','red');
                   
                   xlabel('$x_{1,c}$ [$\mu$m]','fontsize',sizeoffonts);
                   ylabel('$p_{hd,max}$ [MPa]','fontsize',sizeoffonts);   
                   xticks(-500:250:500);
                   yticks(0:250:1500);
                   xlim([-500 500]);
                   ylim([0 1500]);

                   ax = gca;
                   ax.FontSize = sizeofticksfonts;
                   title("$p_{hd, max}$");
                  %==========================================================================
               case 3
                  %========================================================================== 
                  h_surf = subplot(rows, columns, 1);
                   plot_aux.h_ma_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x1;
                   plot_aux.h_ma_surface(it_sim).Nx1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.Nx1;
                   plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x2;
                   plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.Nx2;
                   plot_aux.h_ma_surface(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.h_ma;       

                   %[plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_surface(it_sim).Nx2)/2));       % [m]   x2-position of line plot
                   [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(ceil(plot_aux.h_ma_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

                   surf(plot_aux.h_ma_surface(it_sim).x1*1e6,plot_aux.h_ma_surface(it_sim).x2*1e6,plot_aux.h_ma_surface(it_sim).h_ma'*1e6)
                   hold on
                   plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
                    plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
                    plot_aux.h_ma_surface(it_sim).h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim))'*1e6,'LineWidth',widthlines,'Color',KIT_colorlist{2})

                   xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   zlabel('$h \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   xticks(-500:250:500);
                   yticks(-500:250:500);
                   zlim(ranges.h);
                   fig.CurrentAxes.ZDir = 'Reverse';
                   shading interp
                   light
                   lighting 'gouraud'
                   lightangle(-45,-30)
                   colormap(h_surf,flipud(colmap_h))
                   caxis(ranges.h)
                   ax = gca;
                   ax.FontSize = sizeofticksfonts;
                   title("$h$");

                  %==========================================================================
                  p_surf = subplot(rows, columns, 2);
                   plot_aux.hydr_pressure_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x1;
                   plot_aux.hydr_pressure_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x2;
                   plot_aux.hydr_pressure_surface(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).sol.p_hd;
                   plot_aux.hydr_pressure_surface(it_sim).p_cav = loaded_data(it_sim).fld.p_cav;
                   
                   surf(plot_aux.hydr_pressure_surface(it_sim).x1*1e6,plot_aux.hydr_pressure_surface(it_sim).x2*1e6,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6)

                   hold on
                   contour3(plot_aux.hydr_pressure_surface(it_sim).x1*1e6,plot_aux.hydr_pressure_surface(it_sim).x2*1e6,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6,...
                       [plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6-eps plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6+eps],'LineWidth',widthlines,'Color',KIT_colorlist{3})        
                   hold on
                   plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
                       plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
                       plot_aux.hydr_pressure_surface(it_sim).p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim))'*1e-6,'LineWidth',widthlines,'Color',KIT_colorlist{2})

                   xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   zlabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);
                   xticks(-500:250:500);
                   yticks(-500:250:500);  
                   xlim([-500 500]);
                   ylim([-500 500]);
                   zlim(ranges.p_ax);
                   shading interp
                   light
                   lighting 'gouraud'
                   lightangle(-45,-30)

                   colormap(p_surf,colmap_p)
                   caxis(ranges.p_co)
                   ax = gca;
                   ax.FontSize = sizeofticksfonts;                  
                   title("$p_{hd}$");

                  %==========================================================================
                  load_force = subplot(rows,columns, [3,4]);
                  title("Load Force");

                  %==========================================================================              
               case 4
                  %========================================================================== 
                  h_surf = subplot(rows, columns, 1);
                   plot_aux.h_ma_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x1;
                   plot_aux.h_ma_surface(it_sim).Nx1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.Nx1;
                   plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x2;
                   plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.Nx2;
                   plot_aux.h_ma_surface(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.h_ma;       

                   %[plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_surface(it_sim).Nx2)/2));       % [m]   x2-position of line plot
                   [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(ceil(plot_aux.h_ma_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

                   surf(plot_aux.h_ma_surface(it_sim).x1*1e6,plot_aux.h_ma_surface(it_sim).x2*1e6,plot_aux.h_ma_surface(it_sim).h_ma'*1e6)
                   hold on
                   plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
                    plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
                    plot_aux.h_ma_surface(it_sim).h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim))'*1e6,'LineWidth',widthlines,'Color',KIT_colorlist{2})

                   xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   zlabel('$h \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   xticks(-500:250:500);
                   yticks(-500:250:500);
                   zlim(ranges.h);
                   fig.CurrentAxes.ZDir = 'Reverse';
                   shading interp
                   light
                   lighting 'gouraud'
                   lightangle(-45,-30)
                   colormap(h_surf,flipud(colmap_h))
                   caxis(ranges.h)
                   ax = gca;
                   ax.FontSize = sizeofticksfonts;
                   title("$h$");

                  %==========================================================================
                  p_surf = subplot(rows, columns, 2);
                   plot_aux.hydr_pressure_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x1;
                   plot_aux.hydr_pressure_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x2;
                   plot_aux.hydr_pressure_surface(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).sol.p_hd;
                   plot_aux.hydr_pressure_surface(it_sim).p_cav = loaded_data(it_sim).fld.p_cav;
                   
                   surf(plot_aux.hydr_pressure_surface(it_sim).x1*1e6,plot_aux.hydr_pressure_surface(it_sim).x2*1e6,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6)

                   hold on
                   contour3(plot_aux.hydr_pressure_surface(it_sim).x1*1e6,plot_aux.hydr_pressure_surface(it_sim).x2*1e6,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6,...
                       [plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6-eps plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6+eps],'LineWidth',widthlines,'Color',KIT_colorlist{3})        
                   hold on
                   plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
                       plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
                       plot_aux.hydr_pressure_surface(it_sim).p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim))'*1e-6,'LineWidth',widthlines,'Color',KIT_colorlist{2})

                   xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                   zlabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);
                   xticks(-500:250:500);
                   yticks(-500:250:500);  
                   xlim([-500 500]);
                   ylim([-500 500]);
                   zlim(ranges.p_ax);
                   shading interp
                   light
                   lighting 'gouraud'
                   lightangle(-45,-30)

                   colormap(p_surf,colmap_p)
                   caxis(ranges.p_co)
                   ax = gca;
                   ax.FontSize = sizeofticksfonts;                  
                   title("$p_{hd}$");

                  %==========================================================================               
               case 5
                  %==========================================================================
                  load_force = subplot(rows,columns, [1,2]);
                  title("Load Force");

                  %==========================================================================                    
               case 6
                   if flag_vid.h_surf
                       %========================================================================== 
                       h_surf = subplot(rows, columns, 1);
                        plot_aux.h_ma_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x1;
                        plot_aux.h_ma_surface(it_sim).Nx1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.Nx1;
                        plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x2;
                        plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.Nx2;
                        plot_aux.h_ma_surface(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.h_ma;       

                        %[plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_surface(it_sim).Nx2)/2));       % [m]   x2-position of line plot
                        [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(ceil(plot_aux.h_ma_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

                        surf(plot_aux.h_ma_surface(it_sim).x1*1e6,plot_aux.h_ma_surface(it_sim).x2*1e6,plot_aux.h_ma_surface(it_sim).h_ma'*1e6)
                        hold on
                        plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
                         plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
                         plot_aux.h_ma_surface(it_sim).h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim))'*1e6,'LineWidth',widthlines,'Color',KIT_colorlist{2})

                        xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                        ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                        zlabel('$h \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                        xticks(-500:250:500);
                        yticks(-500:250:500);
                        zlim(ranges.h);
                        fig.CurrentAxes.ZDir = 'Reverse';
                        shading interp
                        light
                        lighting 'gouraud'
                        lightangle(-45,-30)
                        colormap(h_surf,flipud(colmap_h))
                        caxis(ranges.h)
                        ax = gca;
                        ax.FontSize = sizeofticksfonts;
                        title("$h$");

                       %==========================================================================
                   elseif flag_vid.p_surf
                       %========================================================================== 
                       p_surf = subplot(rows, columns, 1);
                        plot_aux.hydr_pressure_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x1;
                        plot_aux.hydr_pressure_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x2;
                        plot_aux.hydr_pressure_surface(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).sol.p_hd;
                        plot_aux.hydr_pressure_surface(it_sim).p_cav = loaded_data(it_sim).fld.p_cav;

                        surf(plot_aux.hydr_pressure_surface(it_sim).x1*1e6,plot_aux.hydr_pressure_surface(it_sim).x2*1e6,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6)

                        hold on
                        contour3(plot_aux.hydr_pressure_surface(it_sim).x1*1e6,plot_aux.hydr_pressure_surface(it_sim).x2*1e6,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6,...
                            [plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6-eps plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6+eps],'LineWidth',widthlines,'Color',KIT_colorlist{3})        
                        hold on
                        plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
                            plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
                            plot_aux.hydr_pressure_surface(it_sim).p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim))'*1e-6,'LineWidth',widthlines,'Color',KIT_colorlist{2})

                        xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                        ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                        zlabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);
                        xticks(-500:250:500);
                        yticks(-500:250:500);  
                        xlim([-500 500]);
                        ylim([-500 500]);
                        zlim(ranges.p_ax);
                        shading interp
                        light
                        lighting 'gouraud'
                        lightangle(-45,-30)

                        colormap(p_surf,colmap_p)
                        caxis(ranges.p_co)
                        ax = gca;
                        ax.FontSize = sizeofticksfonts;                  
                        title("$p_{hd}$");

                       %==========================================================================                   
                   end
               case 7
                   if flag_vid.h_surf
                       %========================================================================== 
                       h_surf = subplot(rows, columns, 1);
                        plot_aux.h_ma_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x1;
                        plot_aux.h_ma_surface(it_sim).Nx1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.Nx1;
                        plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x2;
                        plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.Nx2;
                        plot_aux.h_ma_surface(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.h_ma;       

                        %[plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_surface(it_sim).Nx2)/2));       % [m]   x2-position of line plot
                        [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(ceil(plot_aux.h_ma_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

                        surf(plot_aux.h_ma_surface(it_sim).x1*1e6,plot_aux.h_ma_surface(it_sim).x2*1e6,plot_aux.h_ma_surface(it_sim).h_ma'*1e6)
                        hold on
                        plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
                        plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
                        plot_aux.h_ma_surface(it_sim).h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim))'*1e6,'LineWidth',widthlines,'Color',KIT_colorlist{2})

                        xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                        ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                        zlabel('$h \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                        xticks(-500:250:500);
                        yticks(-500:250:500);
                        zlim(ranges.h);
                        fig.CurrentAxes.ZDir = 'Reverse';
                        shading interp
                        light
                        lighting 'gouraud'
                        lightangle(-45,-30)
                        colormap(h_surf,flipud(colmap_h))
                        caxis(ranges.h)
                        ax = gca;
                        ax.FontSize = sizeofticksfonts;
                        title("$h$");

                       %==========================================================================
                       h_min = subplot(rows, columns, 2);
                        plot(plot_aux.min_h(it_sim).dimple_center*1e6, plot_aux.min_h(it_sim).value*1e9, 'LineWidth',widthlines,'Color',KIT_colorlist{1} );
                        hold on
                        grid on
                        % Red dot:
                        plot(plot_aux.min_h(it_sim).dimple_center(it_time)*1e6, plot_aux.min_h(it_sim).value(it_time)*1e9,'ro','MarkerFaceColor','red','MarkerSize',2);                 
                        % Vertical line:
                        line([plot_aux.min_h(it_sim).dimple_center(it_time)*1e6 plot_aux.min_h(it_sim).dimple_center(it_time)*1e6],[0 plot_aux.min_h(it_sim).value(it_time)*1e9]...
                            ,'LineStyle',':','LineWidth',widthlines/2,'Color','red');                   
                        xlabel('$x_{1,c}$ [$\mu$m]','fontsize',sizeoffonts);
                        ylabel('$h_{ma,min}$ [nm]','fontsize',sizeoffonts);   
                        xticks(-500:250:500);
                        xlim([-500 500]);
                        ylim([0 300]);
 
                        ax = gca;
                        ax.FontSize = sizeofticksfonts;
                        title("$h_{min}$");
                       %==========================================================================
                   elseif flag_vid.p_surf
                       %========================================================================== 
                       p_surf = subplot(rows, columns, 1);
                        plot_aux.hydr_pressure_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x1;
                        plot_aux.hydr_pressure_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.x2;
                        plot_aux.hydr_pressure_surface(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).sol.p_hd;
                        plot_aux.hydr_pressure_surface(it_sim).p_cav = loaded_data(it_sim).fld.p_cav;

                        surf(plot_aux.hydr_pressure_surface(it_sim).x1*1e6,plot_aux.hydr_pressure_surface(it_sim).x2*1e6,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6)

                        hold on
                        contour3(plot_aux.hydr_pressure_surface(it_sim).x1*1e6,plot_aux.hydr_pressure_surface(it_sim).x2*1e6,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6,...
                            [plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6-eps plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6+eps],'LineWidth',widthlines,'Color',KIT_colorlist{3})        
                        hold on
                        plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
                            plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
                            plot_aux.hydr_pressure_surface(it_sim).p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim))'*1e-6,'LineWidth',widthlines,'Color',KIT_colorlist{2})

                        xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                        ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
                        zlabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);
                        xticks(-500:250:500);
                        yticks(-500:250:500);  
                        xlim([-500 500]);
                        ylim([-500 500]);
                        zlim(ranges.p_ax);
                        shading interp
                        light
                        lighting 'gouraud'
                        lightangle(-45,-30)

                        colormap(p_surf,colmap_p)
                        caxis(ranges.p_co)
                        ax = gca;
                        ax.FontSize = sizeofticksfonts;                  
                        title("$p_{hd}$");


                       %==========================================================================
                       p_max = subplot(rows, columns, 2);
                        plot(plot_aux.max_p(it_sim).dimple_center*1e6, plot_aux.max_p(it_sim).value*1e-6, 'LineWidth',widthlines,'Color',KIT_colorlist{1} );
                        hold on
                        grid on
                        % Red dot:
                        plot(plot_aux.max_p(it_sim).dimple_center(it_time)*1e6, plot_aux.max_p(it_sim).value(it_time)*1e-6,'ro','MarkerFaceColor','red','MarkerSize',2);
                        % Vertical line:
                        line([plot_aux.max_p(it_sim).dimple_center(it_time)*1e6 plot_aux.max_p(it_sim).dimple_center(it_time)*1e6],[0 plot_aux.max_p(it_sim).value(it_time)*1e-6]...
                            ,'LineStyle',':','LineWidth',widthlines/2,'Color','red');

                        xlabel('$x_{1,c}$ [$\mu$m]','fontsize',sizeoffonts);
                        ylabel('$p_{hd,max}$ [MPa]','fontsize',sizeoffonts);   
                        xticks(-500:250:500);
                        yticks(0:250:1500);
                        xlim([-500 500]);
                        ylim([0 1500]);
 
                        ax = gca;
                        ax.FontSize = sizeofticksfonts;
                        title("$p_{hd, max}$");
                       %==========================================================================                   
                   end
           end
           sgtitle(loaded_data(it_sim).descr + " @ $t = $ " + loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).alg.time.t*1000 + " $\mathrm{ms}$",'fontsize',sizeoffonts) 
           % Ensure that all frames have the same size by giving enough time until window is reset:
           pause(0.5)
           % Save the frame into a movie vector
           movie_vector(it_time) = getframe(fig); 
       end
       % Save the movie, MPEG-4 is for generating an mp4:
       myWriter = VideoWriter(output_main_path + sprintf("/animation_sim_0%i_config_0%i",it_sim,config),'Uncompressed AVI');
       myWriter.FrameRate = fps;
      
       % Open VideoWriter object, write the movie and close the file:
       open(myWriter);
       writeVideo(myWriter,movie_vector);
       close(myWriter);   
       clear movie_vector;

       
   end
    
end
 fprintf("rows: " + rows + ", columns: " + columns + "\n");
