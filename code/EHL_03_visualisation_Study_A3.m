close all; clc; clearvars; 
% Visualization of the EHL solver results
%
% Recycling of script EHL_03_visualisation_Study_A for visualization of 
% results of simulations with less strict residual definition
%
% Erik Hansen, 04.03.2022
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Execution of this script takes quite long because a lot of data has to be
% read in!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%==========================================================================
%% User input:
%==========================================================================

% Output path:
output_main_path ='./../data/EHL_03_visualisation/Study_A3/';

% Specify if plots should be saved and the plot resolution:
flag_save_plots     = true;                                    % [-]   boolean whether to save the plots and table or not
if flag_save_plots
    print_res       = '-r600';                                  % [-]   resolution used when the plots are printed
end

% Create main output directory:
if flag_save_plots 
    warning('off','MATLAB:MKDIR:DirectoryExists')
    mkdir (output_main_path)
end

% Supply the amount of total simulations:
sim_N                   = 8;

% Supply simulation identifiers:
loaded_data(1).sim_id   = 'Study_A3/Nr_1/';
loaded_data(2).sim_id   = 'Study_A3/Nr_2/';
loaded_data(3).sim_id   = 'Study_A3/Nr_3/';
loaded_data(4).sim_id   = 'Study_A3/Nr_4/';
loaded_data(5).sim_id   = 'Study_A3/Nr_5/';
loaded_data(6).sim_id   = 'Study_A3/Nr_6/';
loaded_data(7).sim_id   = 'Study_A3/Nr_7/';
loaded_data(8).sim_id   = 'Study_A3/Nr_8/';

% Supply simulation descriptions:
loaded_data(1).descr    = 'Deep, UI, SSR = 0';
loaded_data(2).descr    = 'Shallow, UI, SSR = 0';
loaded_data(3).descr    = 'Deep, UI, SSR = -0.5';
loaded_data(4).descr    = 'Shallow, UI, SSR = -0.5';
loaded_data(5).descr    = 'Deep, QUICK, SSR = 0';
loaded_data(6).descr    = 'Shallow, QUICK, SSR = 0';
loaded_data(7).descr    = 'Deep, QUICK, SSR = -0.5';
loaded_data(8).descr    = 'Shallow, QUICK, SSR = -0.5';

% Choose which simulations should be plotted:
sims_to_plot = [ 1 2 3 4 5 6 7 8];

% Choose detailed operating condition:
i_OC          = 1;

% Choose time point to plot in detail:
i_time        = 1;

% Define flags for creating the figures in the paper:
flag.mourier_comp = true; % generates a comparison figure containing the results from Mourier et al. (2006) 
flag.cav_comp     = true;
flag.p_comp       = true;
flag.h_comp       = true;
flag.titles_on    = false;

% Define general flags for creating plots:
flag.h_surf       = false;
flag.h_line       = false;
flag.p_surf       = false;
flag.p_line       = false;
flag.FBNS_res     = false;
flag.load_bal_res = false;
flag.h_min        = false;
flag.p_max        = false;
 

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
% Loading in data from Mourier 2006:
%==========================================================================
% The experimental results only have data about the gap height, the
% simulation results have the pressure distribution additionally:
if flag.mourier_comp
    % Initialize the cell arrays to hold the resulting plots:
    % For SSR = 0:
    Mourier_exp = cell(5,1);
    Mourier_sim = cell(5,2); % First column is the gap height, the second the pressure. The rows correspond to different time points
    % For SSR = -0.5:
    Mourier_2_exp  = cell(5,1);
    Mourier_2_sim  = cell(5,2);
    % Read in the data:
    % For SSR=0:
    Mourier_exp{1,1} = readmatrix('../data/Mourier_2006/mourier_ssr0_exp_1.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_exp{2,1} = readmatrix('../data/Mourier_2006/mourier_ssr0_exp_2.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_exp{3,1} = readmatrix('../data/Mourier_2006/mourier_ssr0_exp_3.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_exp{4,1} = readmatrix('../data/Mourier_2006/mourier_ssr0_exp_4.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_exp{5,1} = readmatrix('../data/Mourier_2006/mourier_ssr0_exp_5.csv','NumHeaderLines',1,'DecimalSeparator','.');

    Mourier_sim{1,1} =  readmatrix('../data/Mourier_2006/mourier_ssr0_sim_1_h.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_sim{2,1} =  readmatrix('../data/Mourier_2006/mourier_ssr0_sim_2_h.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_sim{3,1} =  readmatrix('../data/Mourier_2006/mourier_ssr0_sim_3_h.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_sim{4,1} =  readmatrix('../data/Mourier_2006/mourier_ssr0_sim_4_h.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_sim{5,1} =  readmatrix('../data/Mourier_2006/mourier_ssr0_sim_5_h.csv','NumHeaderLines',1,'DecimalSeparator','.');

    Mourier_sim{1,2} =  readmatrix('../data/Mourier_2006/mourier_ssr0_sim_1_p.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_sim{2,2} =  readmatrix('../data/Mourier_2006/mourier_ssr0_sim_2_p.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_sim{3,2} =  readmatrix('../data/Mourier_2006/mourier_ssr0_sim_3_p.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_sim{4,2} =  readmatrix('../data/Mourier_2006/mourier_ssr0_sim_4_p.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_sim{5,2} =  readmatrix('../data/Mourier_2006/mourier_ssr0_sim_5_p.csv','NumHeaderLines',1,'DecimalSeparator','.');    
    
    % For SSR=-0.5:
    Mourier_2_exp{1,1} = readmatrix('../data/Mourier_2006/mourier_ssr-05_exp_1.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_2_exp{2,1} = readmatrix('../data/Mourier_2006/mourier_ssr-05_exp_2.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_2_exp{3,1} = readmatrix('../data/Mourier_2006/mourier_ssr-05_exp_3.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_2_exp{4,1} = readmatrix('../data/Mourier_2006/mourier_ssr-05_exp_4.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_2_exp{5,1} = readmatrix('../data/Mourier_2006/mourier_ssr-05_exp_5.csv','NumHeaderLines',1,'DecimalSeparator','.');    
    
    Mourier_2_sim{1,1} =  readmatrix('../data/Mourier_2006/mourier_ssr-05_sim_1_h.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_2_sim{2,1} =  readmatrix('../data/Mourier_2006/mourier_ssr-05_sim_2_h.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_2_sim{3,1} =  readmatrix('../data/Mourier_2006/mourier_ssr-05_sim_3_h.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_2_sim{4,1} =  readmatrix('../data/Mourier_2006/mourier_ssr-05_sim_4_h.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_2_sim{5,1} =  readmatrix('../data/Mourier_2006/mourier_ssr-05_sim_5_h.csv','NumHeaderLines',1,'DecimalSeparator','.');

    Mourier_2_sim{1,2} =  readmatrix('../data/Mourier_2006/mourier_ssr-05_sim_1_p.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_2_sim{2,2} =  readmatrix('../data/Mourier_2006/mourier_ssr-05_sim_2_p.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_2_sim{3,2} =  readmatrix('../data/Mourier_2006/mourier_ssr-05_sim_3_p.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_2_sim{4,2} =  readmatrix('../data/Mourier_2006/mourier_ssr-05_sim_4_p.csv','NumHeaderLines',1,'DecimalSeparator','.');
    Mourier_2_sim{5,2} =  readmatrix('../data/Mourier_2006/mourier_ssr-05_sim_5_p.csv','NumHeaderLines',1,'DecimalSeparator','.');    
    
    % Convert units to [m] and [Pa]:
    for outer_row = 1:5
        % For SSR = :0
        Mourier_exp{outer_row,1}(:,1) = Mourier_exp{outer_row,1}(:,1)*1e-6; % converting the spatial coordinate from micrometers to meters
        Mourier_exp{outer_row,1}(:,2) = Mourier_exp{outer_row,1}(:,2)*1e-9; % converting the gap height from nanometers to meters
        
        Mourier_sim{outer_row,1}(:,1) = Mourier_sim{outer_row,1}(:,1)*1e-6; % converting the spatial coordinate from micrometers to meters
        Mourier_sim{outer_row,1}(:,2) = Mourier_sim{outer_row,1}(:,2)*1e-9; % converting the gap height from nanometers to meters

        Mourier_sim{outer_row,2}(:,1) = Mourier_sim{outer_row,2}(:,1)*1e-6; % converting the spatial coordinate from micrometers to meters for the pressure plot
        Mourier_sim{outer_row,2}(:,2) = Mourier_sim{outer_row,2}(:,2)*1e6; % converting the pressure from MPa to Pa   
        
        % For SSR=-0.5:
        Mourier_2_exp{outer_row,1}(:,1) = Mourier_2_exp{outer_row,1}(:,1)*1e-6; % converting the spatial coordinate from micrometers to meters
        Mourier_2_exp{outer_row,1}(:,2) = Mourier_2_exp{outer_row,1}(:,2)*1e-9; % converting the gap height from nanometers to meters
        
        Mourier_2_sim{outer_row,1}(:,1) = Mourier_2_sim{outer_row,1}(:,1)*1e-6; % converting the spatial coordinate from micrometers to meters
        Mourier_2_sim{outer_row,1}(:,2) = Mourier_2_sim{outer_row,1}(:,2)*1e-9; % converting the gap height from nanometers to meters

        Mourier_2_sim{outer_row,2}(:,1) = Mourier_2_sim{outer_row,2}(:,1)*1e-6; % converting the spatial coordinate from micrometers to meters for the pressure plot
        Mourier_2_sim{outer_row,2}(:,2) = Mourier_2_sim{outer_row,2}(:,2)*1e6; % converting the pressure from MPa to Pa              
    end   
end

% Create a "analysis" object to reference dimple positions and times with
% time point indices:
% Initialize the object:
analysis(1).sim_id   = 'Study_A/Nr_1/';
analysis(2).sim_id   = 'Study_A/Nr_2/';
analysis(3).sim_id   = 'Study_A/Nr_3/';
analysis(4).sim_id   = 'Study_A/Nr_4/';
analysis(5).sim_id   = 'Study_A/Nr_5/';
analysis(6).sim_id   = 'Study_A/Nr_6/';
analysis(7).sim_id   = 'Study_A/Nr_7/';
analysis(8).sim_id   = 'Study_A/Nr_8/';

for it_sim = 1:sim_N
    if ismember(it_sim,sims_to_plot)
        % Get how many time steps are used in this simulation:
        time_N = loaded_data(it_sim).loaded_OC.loaded_Time(1).alg.time.N;
        % Initialize the reference matrix:
        analysis(it_sim).reference = zeros(time_N + 1,3); % columns: time index, time, dimple x1 position
        % Name the columns:
        analysis(it_sim).reference(1,1) = "time index"; 
        analysis(it_sim).reference(1,2) = "time"; 
        analysis(it_sim).reference(1,3) = "dimple x1 position";           
        for it_time = 1:time_N
            % Current time:
            current_time = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).alg.time.t;  
            % Current dimple position:
            current_position = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).slc.x1_c;  
            % Save the values:
            analysis(it_sim).reference(it_time + 1,1) = it_time; 
            analysis(it_sim).reference(it_time + 1,2) = current_time; 
            analysis(it_sim).reference(it_time + 1,3) = current_position;             
        end
    end 
end
 clear current_position; clear current_time; clear time_N;
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
colmap_h = batlow;
colmap_p = batlow;

batlowS_list = cell(1,100);
for i = 1:100
    batlowS_list{i} = [batlowS(i,1), batlowS(i,2), batlowS(i,3)];
end

KIT_colorlist = batlowS_list(3:10);

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
ranges.h_flag = false;
if ranges.h_flag
    ranges.h = [0 5e-6]; % [m]
end
ranges.x1_flag = false;
if ranges.x1_flag
    ranges.x1 = [-6e-3 6e-3]; % [m]
end
ranges.x2_flag = false;
if ranges.x2_flag
    ranges.x2 = [-6e-3 6e-3]; % [m]
end
ranges.it_tot_flag = false;
if ranges.it_tot_flag
   ranges.it_tot = [0 500]; % [-] 
end
ranges.min_h_flag = false;
if ranges.min_h_flag
   ranges.min_h = [0 100e-9]; % [m]
end
ranges.max_p_flag = false;
if ranges.max_p_flag
   ranges.max_p = [800e6 1e9]; % [Pa]
end


%% Plotting - Macroscopic gap height:
%==========================================================================
% Surface:
%==========================================================================
if flag.h_surf
    for it_sim = 1:sim_N
       if ismember(it_sim,sims_to_plot)
           fig = figure('Units','centimeters','Position',[1 1 figuresize]);
           plot_aux.h_ma_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
           plot_aux.h_ma_surface(it_sim).Nx1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.Nx1;
           plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x2;
           plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.Nx2;
           plot_aux.h_ma_surface(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.h_ma;       

           %[plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_surface(it_sim).Nx2)/2));       % [m]   x2-position of line plot
           [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(ceil(plot_aux.h_ma_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

           surf(plot_aux.h_ma_surface(it_sim).x1*1e6,plot_aux.h_ma_surface(it_sim).x2*1e6,plot_aux.h_ma_surface(it_sim).h_ma'*1e9)
           hold on
           plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
            plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
            plot_aux.h_ma_surface(it_sim).h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim))'*1e9,'LineWidth',widthlines,'Color',KIT_colorlist{1})

           xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
           ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
           zlabel('$h_{ma} \mathrm{[nm]}$','fontsize',sizeoffonts);

           fig.CurrentAxes.ZDir = 'Reverse';
           shading interp
           light
           lighting 'gouraud'
           lightangle(-45,-30)
           colormap(flipud(colmap_h))


           ax = gca;
           ax.FontSize = sizeofticksfonts;
           if flag.titles_on           
               title(loaded_data(it_sim).descr + " gap height @ $U =$ " + loaded_data(it_sim).opc.u_up(i_OC) + "$\mathrm{m/s}$",'fontsize',sizeoffonts)
           end
           if ranges.x1_flag
               xlim(ranges.x1*1e6)
           end
           if ranges.x2_flag
               ylim(ranges.x2*1e6)
           end
           if ranges.h_flag
               caxis(ranges.h*1e9)
               zlim(ranges.h*1e9)
           end
           if flag_save_plots        
               plotname = char(compose('h_ma_surface_%i_%i',it_sim,i_OC));
               print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
               print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
               savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
               clear plotname;
           end
           clear fig; clear ax;
       end
    end
end
%==========================================================================
% Line plots:
%==========================================================================
if flag.h_line
    for it_sim = 1:sim_N 
       if ismember(it_sim, sims_to_plot) 
        plot_aux.h_ma_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
        plot_aux.h_ma_line(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim));

        fig = figure('Units','centimeters','Position',[21.1 1 figuresize]);
        plot(plot_aux.h_ma_line(it_sim).x1*1e6, ...
            plot_aux.h_ma_line(it_sim).h_ma*1e9,...
            'LineWidth',widthlines,'Color',KIT_colorlist{1})
        grid on
        fig.CurrentAxes.YDir = 'Reverse';
        xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
        ylabel('$h_{ma} \mathrm{[nm]}$','fontsize',sizeoffonts);
        if ranges.x1_flag
            xlim(ranges.x1*1e6)
        end
        if ranges.h_flag
            ylim(ranges.h*1e9)
        end
        ax = gca;
        ax.FontSize = sizeofticksfonts;
        if flag.titles_on        
            title(loaded_data(it_sim).descr + " gap height line plot @$x_2= $" + plot_aux.h_ma_line_x2_0(it_sim)*1e3 + "$\mathrm{mm}$",'fontsize',sizeoffonts)
        end
        if flag_save_plots        
            plotname = char(compose('h_ma_line_%i_%i',it_sim,i_OC));
            print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
            print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
            savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
            clear plotname;
        end
        clear fig; clear ax;
       end
    end
end
%% Plotting - Hydrodynamic pressure:
%==========================================================================
% Surface:
%==========================================================================
if flag.p_surf
    for it_sim = 1:sim_N 
       if ismember(it_sim, sims_to_plot) 
        fig = figure('Units','centimeters','Position',[11.1 11 figuresize]);

        plot_aux.hydr_pressure_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
        plot_aux.hydr_pressure_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x2;
        plot_aux.hydr_pressure_surface(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd;
        plot_aux.hydr_pressure_surface(it_sim).p_cav = loaded_data(it_sim).fld.p_cav;

        surf(plot_aux.hydr_pressure_surface(it_sim).x1*1e6,plot_aux.hydr_pressure_surface(it_sim).x2*1e6,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6)

        hold on
        contour3(plot_aux.hydr_pressure_surface(it_sim).x1*1e6,plot_aux.hydr_pressure_surface(it_sim).x2*1e6,plot_aux.hydr_pressure_surface(it_sim).p_hd'*1e-6,...
            [plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6-eps plot_aux.hydr_pressure_surface(it_sim).p_cav*1e-6+eps],'LineWidth',widthlines,'Color',KIT_colorlist{2})        
        hold on
        plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
            plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
            plot_aux.hydr_pressure_surface(it_sim).p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim))'*1e-6,'LineWidth',widthlines,'Color',KIT_colorlist{1})

        xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
        ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
        zlabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);

        shading interp
        light
        lighting 'gouraud'
        lightangle(-45,-30)

        colormap(colmap_p)
        ax = gca;
        ax.FontSize = sizeofticksfonts;
        if flag.titles_on        
            title(loaded_data(it_sim).descr + " pressure @ $U =$ " + loaded_data(it_sim).opc.u_up(i_OC) + "$\mathrm{m/s}$",'fontsize',sizeoffonts)
        end
        if ranges.x1_flag
            xlim(ranges.x1*1e6)
        end
        if ranges.x2_flag
            ylim(ranges.x2*1e6)
        end
        if ranges.p_flag
            caxis(ranges.p*1e-6)
            zlim(ranges.p*1e-6)
        end
        if flag_save_plots        
            plotname = char(compose('hydr_pressure_surface_%i_%i',it_sim,i_OC));
            print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
            print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
            savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
            clear plotname;
        end
        clear fig; clear ax;
       end
    end
end
%==========================================================================
% Line plots:
%==========================================================================
if flag.p_line
    for it_sim = 1:sim_N 
       if ismember(it_sim, sims_to_plot) 
        plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
        plot_aux.hydr_pressure_line(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim));

        fig = figure('Units','centimeters','Position',[11.1 1 figuresize]);
        plot(plot_aux.hydr_pressure_line(it_sim).x1*1e6, ...
            plot_aux.hydr_pressure_line(it_sim).p_hd*1e-6,...
            'LineWidth',widthlines,'Color',KIT_colorlist{1})
        grid on
        xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
        ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);
        if ranges.x1_flag
            xlim(ranges.x1*1e6)
        end
        if ranges.h_flag
            ylim(ranges.p*1e6)
        end
        ax = gca;
        ax.FontSize = sizeofticksfonts;
        if flag.titles_on        
            title(loaded_data(it_sim).descr + " pressure line plot @$x_2= $" + plot_aux.h_ma_line_x2_0(it_sim)*1e3 + "$\mathrm{mm}$",'fontsize',sizeoffonts)
        end
        if flag_save_plots
            plotname = char(compose('hydr_pressure_line_%i_%i',it_sim,i_OC));
            print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
            print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
            savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
            clear plotname;
        end
        clear fig; clear ax;    
       end
    end
end
%% Plotting - Residuals:
%==========================================================================
% FBNS Residuals:
%==========================================================================
if flag.FBNS_res
    for it_sim = 1:sim_N
       if ismember(it_sim,sims_to_plot)
           fig = figure('Units','centimeters','Position',[15.1 11 figuresize]);
           iterations     = 1:loaded_data(it_sim).loaded_OC.loaded_Time(i_time).alg.it_tot;
           FBNS_residuals = loaded_data(it_sim).loaded_OC.loaded_Time(i_time).res.FBNS.FBNS;
           semilogy(iterations, FBNS_residuals,'LineWidth',widthlines,'Color',KIT_colorlist{1});

           grid on
           xlabel("Total iterations for this time point. [-]",'fontsize',sizeoffonts);
           ylabel("Chosen FBNS residual. [-]",'fontsize',sizeoffonts)
           if ranges.it_tot_flag
              xlim(ranges.it_tot); 
           end
           ax = gca;
           ax.FontSize = sizeofticksfonts;
           if flag.titles_on                          
               title(loaded_data(it_sim).descr + " relative FBNS residual",'fontsize',sizeoffonts);
           end
           if flag_save_plots        
               plotname = char(compose('FBNS_res_%i_%i',it_sim,i_OC));
               print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
               print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
               savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
               clear plotname;
           end
        clear fig; clear ax; clear FBNS_residuals; clear iterations;
       end
    end
end 
%==========================================================================
% Load Residuals:
%==========================================================================
if flag.load_bal_res
    for it_sim = 1:sim_N
       if ismember(it_sim,sims_to_plot)
           fig = figure('Units','centimeters','Position',[15.1 1 figuresize]);
           iterations     = 1:loaded_data(it_sim).loaded_OC.loaded_Time(i_time).alg.it_tot;
           Load_residuals = loaded_data(it_sim).loaded_OC.loaded_Time(i_time).res.load.F_N;
           semilogy(iterations, Load_residuals,'LineWidth',widthlines,'Color',KIT_colorlist{1});

           grid on
           warning('off','MATLAB:Axes:NegativeDataInLogAxis'); % NOTE TO ERIK: I think the FBNS residual gets so small that MATLAB thinks the values are negative and gives a warning. I supressed that, it says it just ignores the negative data. :)
           xlabel("Total iterations for this time point. [-]",'fontsize',sizeoffonts);
           ylabel("Chosen Load residual. [-]",'fontsize',sizeoffonts)
           if ranges.it_tot_flag
              xlim(ranges.it_tot); 
           end
           ax = gca;
           ax.FontSize = sizeofticksfonts;
           if flag.titles_on           
               title(loaded_data(it_sim).descr + " relative Load residual",'fontsize',sizeoffonts);
           end
           if flag_save_plots        
               plotname = char(compose('Load_res_%i_%i',it_sim,i_OC));
               print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
                   print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
                   savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
                   clear plotname;
               end
            clear fig; clear ax; clear Load_residuals; clear iterations;
       end
    end
end
%% Plotting - Min Gap height, Max Pressure vs. Dimple Center Position
%==========================================================================
% Minimum gap height:
%==========================================================================
if flag.h_min
    for it_sim = 1:sim_N
       if ismember(it_sim,sims_to_plot)
           fig = figure('Units','centimeters','Position',[5.1 11 figuresize]);
           % Wrapper to go through all time points for this OC:
           total_time_points = loaded_data(it_sim).opc.N_t;
           for it_time = 1:total_time_points
               h_for_this_time    = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).h.h_ma;
               x1_c_for_this_time = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).slc.x1_c;  
               plot_aux.min_h(it_sim).value(it_time)         = min(min(h_for_this_time));
               plot_aux.min_h(it_sim).dimple_center(it_time) = x1_c_for_this_time;          
           end

           plot(plot_aux.min_h(it_sim).dimple_center*1e6, plot_aux.min_h(it_sim).value*1e9, 'LineWidth',widthlines,'Color',KIT_colorlist{1} );
           xlabel('$x_1$-Coordinate of the micro-cavity center [$\mu$m]','fontsize',sizeoffonts);
           ylabel('Minimum gap height ($h_{ma,min}$) across the contact [nm]','fontsize',sizeoffonts);

           if ranges.min_h_flag
              ylim(ranges.min_h*1e9);
           end
           grid on
           ax = gca;
           ax.FontSize = sizeofticksfonts;
           if flag.titles_on           
               title(loaded_data(it_sim).descr + " minimum gap height",'fontsize',sizeoffonts);   
           end
           if flag_save_plots        
               plotname = char(compose('h_min_%i_%i',it_sim,i_OC));
               print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) ;
               print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res);
               savefig(fig,fullfile(output_main_path,append(plotname,'.fig')));
               clear plotname;
           end
           clear fig; clear ax; clear total_time_points; clear x1_c_for_this_time; clear h_for_this_time;
       end
    end
end
%==========================================================================
% Maximum pressure:
%==========================================================================
if flag.p_max
    for it_sim = 1:sim_N
       if ismember(it_sim,sims_to_plot)
           fig = figure('Units','centimeters','Position',[5.1 1 figuresize]);
           % Wrapper to go through all time points for this OC:
           total_time_points = loaded_data(it_sim).opc.N_t;
           for it_time = 1:total_time_points
               p_for_this_time    = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).sol.p_hd  ;
               x1_c_for_this_time = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).slc.x1_c;  
               plot_aux.max_p(it_sim).value(it_time)         = max(max(p_for_this_time));
               plot_aux.max_p(it_sim).dimple_center(it_time) = x1_c_for_this_time;          
           end

           plot(plot_aux.max_p(it_sim).dimple_center*1e6, plot_aux.max_p(it_sim).value*1e-6, 'LineWidth',widthlines,'Color',KIT_colorlist{1} );
           xlabel('$x_1$-Coordinate of the micro-cavity center [$\mu$m]','fontsize',sizeoffonts);
           ylabel('Maximum hydrodynamic pressure ($p_{hd,max}$) across the contact [MPa]','fontsize',sizeoffonts);

           if ranges.max_p_flag
              ylim(ranges.max_p*1e-6);
           end
           grid on
           ax = gca;
           ax.FontSize = sizeofticksfonts;
           if flag.titles_on           
               title(loaded_data(it_sim).descr + " maximum pressure",'fontsize',sizeoffonts);   
           end
           if flag_save_plots        
               plotname = char(compose('p_max_%i_%i',it_sim,i_OC));
               print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) ;
               print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res);
               savefig(fig,fullfile(output_main_path,append(plotname,'.fig')));
               clear plotname;
           end
           clear fig; clear ax; clear total_time_points; clear p_for_this_time; clear x1_c_for_this_time;
       end
    end
end
%% Paper - Surface plots:
fig = figure('Units','centimeters','Position',[11.1 1 13 figuresize(2)*3]);
%==========================================================================
% Macroscopic gap height for UI:
%==========================================================================
it_sim = 1;
su_pl_1 = subplot(3,2,1);
   plot_aux.h_ma_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
   plot_aux.h_ma_surface(it_sim).Nx1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.Nx1;
   plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x2;
   plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.Nx2;
   plot_aux.h_ma_surface(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.h_ma;       

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
   fig.CurrentAxes.ZDir = 'Reverse';
   shading interp
   light
   lighting 'gouraud'
   lightangle(-45,-30)
   colormap(su_pl_1,flipud(colmap_h))


   ax = gca;
   ax.FontSize = sizeofticksfonts;
   if flag.titles_on           
       title(loaded_data(it_sim).descr + " gap height @ $U =$ " + loaded_data(it_sim).opc.u_up(i_OC) + "$\mathrm{m/s}$",'fontsize',sizeoffonts)
    else
        title("\textbf{(a)}")   
   end
   if ranges.x1_flag
       xlim(ranges.x1*1e6)
   end
   if ranges.x2_flag
       ylim(ranges.x2*1e6)
   end
   if ranges.h_flag
       caxis(ranges.h*1e6)
       zlim(ranges.h*1e6)
   end
%==========================================================================
% Macroscopic gap height for QUICK:
%==========================================================================
it_sim = 5;
su_pl_2 = subplot(3,2,2);
   plot_aux.h_ma_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
   plot_aux.h_ma_surface(it_sim).Nx1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.Nx1;
   plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x2;
   plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.Nx2;
   plot_aux.h_ma_surface(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.h_ma;       

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
   
   fig.CurrentAxes.ZDir = 'Reverse';
   shading interp
   light
   lighting 'gouraud'
   lightangle(-45,-30)
   colormap(su_pl_2,flipud(colmap_h))


   ax = gca;
   ax.FontSize = sizeofticksfonts;
   if flag.titles_on           
       title(loaded_data(it_sim).descr + " gap height @ $U =$ " + loaded_data(it_sim).opc.u_up(i_OC) + "$\mathrm{m/s}$",'fontsize',sizeoffonts)
    else
        title("\textbf{(b)}")   
   end
   if ranges.x1_flag
       xlim(ranges.x1*1e6)
   end
   if ranges.x2_flag
       ylim(ranges.x2*1e6)
   end
   if ranges.h_flag
       caxis(ranges.h*1e6)
       zlim(ranges.h*1e6)
   end   
%==========================================================================
% Hydrodynamic pressure for UI:
%==========================================================================
it_sim = 1;
su_pl_3 = subplot(3,2,3);
    plot_aux.hydr_pressure_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.hydr_pressure_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x2;
    plot_aux.hydr_pressure_surface(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd;
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
    zlim([0 400]);
    shading interp
    light
    lighting 'gouraud'
    lightangle(-45,-30)

    colormap(su_pl_3,colmap_p)
    ax = gca;
    ax.FontSize = sizeofticksfonts;
    if flag.titles_on        
        title(loaded_data(it_sim).descr + " pressure @ $U =$ " + loaded_data(it_sim).opc.u_up(i_OC) + "$\mathrm{m/s}$",'fontsize',sizeoffonts)
    else
        title("\textbf{(c)}")
    end
    if ranges.x1_flag
        xlim(ranges.x1*1e6)
    end
    if ranges.x2_flag
        ylim(ranges.x2*1e6)
    end
    if ranges.p_flag
        caxis(ranges.p*1e-6)
        zlim(ranges.p*1e-6)
    end
%==========================================================================
% Hydrodynamic pressure for QUICK:
%==========================================================================
it_sim = 5;
su_pl_4 = subplot(3,2,4);
    plot_aux.hydr_pressure_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.hydr_pressure_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x2;
    plot_aux.hydr_pressure_surface(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.p_hd;
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
    zlim([0 400]);
    shading interp
    light
    lighting 'gouraud'
    lightangle(-45,-30)

    colormap(su_pl_4,colmap_p)
    ax = gca;
    ax.FontSize = sizeofticksfonts;
    if flag.titles_on        
        title(loaded_data(it_sim).descr + " pressure @ $U =$ " + loaded_data(it_sim).opc.u_up(i_OC) + "$\mathrm{m/s}$",'fontsize',sizeoffonts)
    else
        title("\textbf{(d)}")    
    end
    if ranges.x1_flag
        xlim(ranges.x1*1e6)
    end
    if ranges.x2_flag
        ylim(ranges.x2*1e6)
    end
    if ranges.p_flag
        caxis(ranges.p*1e-6)
        zlim(ranges.p*1e-6)
    end    
%==========================================================================
% Cavity fraction for UI:
%==========================================================================
it_sim = 1;
su_pl_5 = subplot(3,2,5);
    plot_aux.theta_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.theta_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x2;
    plot_aux.theta_surface(it_sim).theta = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.thet;
    
    surf(plot_aux.theta_surface(it_sim).x1*1e6,plot_aux.theta_surface(it_sim).x2*1e6,plot_aux.theta_surface(it_sim).theta')
    
    hold on
    plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
        plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
        plot_aux.theta_surface(it_sim).theta(:,plot_aux.h_ma_line_x2_0_i(it_sim))','LineWidth',widthlines,'Color',KIT_colorlist{2})
    
    xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    zlabel('$\theta \mathrm{[-]}$','fontsize',sizeoffonts);
    xticks(-500:250:500);
    yticks(-500:250:500);    
    
    shading interp
    light
    lighting 'gouraud'
    lightangle(-45,-30)
    
    colormap(su_pl_5,colmap_p)
    ax = gca;
    ax.FontSize = sizeofticksfonts;
    if flag.titles_on    
        title(loaded_data(it_sim).descr + " $\theta$ @ $U =$ " + loaded_data(it_sim).opc.u_up(i_OC) + "$\mathrm{m/s}$",'fontsize',sizeoffonts)
    else
        title("\textbf{(e)}")    
    end
    if ranges.x1_flag
        xlim(ranges.x1*1e6)
    end
    if ranges.x2_flag
        ylim(ranges.x2*1e6)
    end
    if ranges.p_flag
        caxis([0 0.5])
        zlim([0 0.5])
    end

%==========================================================================
% Cavity fraction for QUICK:
%==========================================================================
it_sim = 5;
su_pl_6 = subplot(3,2,6);
    plot_aux.theta_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x1;
    plot_aux.theta_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).h.x2;
    plot_aux.theta_surface(it_sim).theta = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(i_time).sol.thet;
    
    surf(plot_aux.theta_surface(it_sim).x1*1e6,plot_aux.theta_surface(it_sim).x2*1e6,plot_aux.theta_surface(it_sim).theta')
    
    hold on
    plot3(plot_aux.h_ma_surface(it_sim).x1*1e6,...
        plot_aux.h_ma_surface(it_sim).x2(plot_aux.h_ma_line_x2_0_i(it_sim))*ones(plot_aux.h_ma_surface(it_sim).Nx1,1)*1e6,...
        plot_aux.theta_surface(it_sim).theta(:,plot_aux.h_ma_line_x2_0_i(it_sim))','LineWidth',widthlines,'Color',KIT_colorlist{2})
    
    xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    zlabel('$\theta \mathrm{[-]}$','fontsize',sizeoffonts);
    xticks(-500:250:500);
    yticks(-500:250:500);
   
    shading interp
    light
    lighting 'gouraud'
    lightangle(-45,-30)
    
    colormap(su_pl_6,colmap_p)
    ax = gca;
    ax.FontSize = sizeofticksfonts;
    if flag.titles_on    
        title(loaded_data(it_sim).descr + " $\theta$ @ $U =$ " + loaded_data(it_sim).opc.u_up(i_OC) + "$\mathrm{m/s}$",'fontsize',sizeoffonts)
    else
        title("\textbf{(f)}")    
    end
    if ranges.x1_flag
        xlim(ranges.x1*1e6)
    end
    if ranges.x2_flag
        ylim(ranges.x2*1e6)
    end
    if ranges.p_flag
        caxis([0 0.5])
        zlim([0 0.5])
    end
if flag_save_plots        
    plotname = char(compose('paper_surfaces'));
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;
%% Paper - Mourier Comparision SSR = 0 Deep Dimple
%==========================================================================
% Deep dimple:
%==========================================================================
% Choose which time points to plot:
times_to_plot = [73 97 116 135 163];
if flag.mourier_comp
    % Define figure:
    fig = figure('Units','Centimeters','Position', [0.1 0.1 34 14]);
    % Define legend:
    mourier_deep_legend = cell(3,1);
    mourier_deep_legend{1} = "EHL-FBNS UI";
    mourier_deep_legend{2} = "EHL-FBNS QUICK";
    mourier_deep_legend{3} = "Experiment - Mourier et al."; 

    

    for plot_index = 1:length(times_to_plot)
        time_to_plot = times_to_plot(plot_index);
        if plot_index == 2
           time_to_plot = 101; %since the sampling coordinate from Mourier is not precise for the second plot a different time point is needed compared to the simulations 
        end
        subplot('Position',[0.05+(plot_index-1)*(0.15+0.05) 0.60 0.12 0.30]);
        hold on
        % Gap height:
        %==========================================================================
        % 1st Order Upwind:
        %==========================================================================
        it_sim = 1;
        % Gap height:
        plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x2;
        plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.Nx2;
        [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(ceil(plot_aux.h_ma_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

        plot_aux.h_ma_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.h_ma_line(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim));

        h_exp_1 = plot(plot_aux.h_ma_line(it_sim).x1*1e6, ...
                   plot_aux.h_ma_line(it_sim).h_ma*1e9,...
                   'LineWidth',widthlines,'Color',KIT_colorlist{1});
        grid on
        xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
        ylabel('$h \mathrm{[nm]}$','fontsize',sizeoffonts);  
        title(['$x_{1,d} =$ ', num2str(analysis(it_sim).reference(1+time_to_plot,3)*1e6,'%.2f'), ' $\mathrm{\mu m}$'],'fontsize',sizeoffonts)
        
        %==========================================================================
        % 2nd Order Quick:
        %==========================================================================
        it_sim = 5;
        
        plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x2;
        plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.Nx2;
        [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(ceil(plot_aux.h_ma_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

        plot_aux.h_ma_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.h_ma_line(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim));
        
        h_exp_2 = plot(plot_aux.h_ma_line(it_sim).x1*1e6, ...
                   plot_aux.h_ma_line(it_sim).h_ma*1e9,...
                   'LineWidth',widthlines,'Color',KIT_colorlist{2},'LineStyle','--');
        %==========================================================================
        % Mourier et al. (2006):
        %==========================================================================
        
        plot_aux.Mourier_exp.x1 = Mourier_exp{plot_index,1}(:,1);
        plot_aux.Mourier_exp.h_ma = Mourier_exp{plot_index,1}(:,2);
        
        h_exp_3 = plot(plot_aux.Mourier_exp.x1*1e6, plot_aux.Mourier_exp.h_ma*1e9,'LineWidth', widthlines,'Color',KIT_colorlist{3},'LineStyle',':');
        set(gca, 'YDir','reverse'); 
        %==========================================================================
        % Defining ranges:
        %==========================================================================
        xlim([-200 200]); % [micrometer]
        xticks([-200 -100 0 100 200])
        ylim([0 350]); %[nanometer]        
        yticks([0 50 100 150 200 250 300 350])
        hold off
        
        % Pressures:
        subplot('Position',[0.05+(plot_index-1)*(0.15+0.05) 0.20 0.12 0.30]);
        hold on
        grid on
        %==========================================================================
        % 1st Order Upwind:
        %==========================================================================        
        it_sim = 1;
        plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.hydr_pressure_line(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).sol.p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim));
        
        p_exp_1 = plot(plot_aux.hydr_pressure_line(it_sim).x1*1e6, ...
              plot_aux.hydr_pressure_line(it_sim).p_hd*1e-6,...
              'LineWidth',widthlines,'Color',KIT_colorlist{1});
        xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
        ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);     
        
        %==========================================================================
        % 2. Order QUICK
        %==========================================================================
        it_sim = 5;
        plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.hydr_pressure_line(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).sol.p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim));
        
        p_exp_2 = plot(plot_aux.hydr_pressure_line(it_sim).x1*1e6, ...
              plot_aux.hydr_pressure_line(it_sim).p_hd*1e-6,...
              'LineWidth',widthlines,'Color',KIT_colorlist{2},'LineStyle','--');
        ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);  
        % Defining ranges:
        xlim([-200 200]); % [micrometer]  
        xticks([-200 -100 0 100 200])
        ylim([0 800]); %[MPa]
        hold off
    end
    % Creating the legend:
    leg = legend([h_exp_1, h_exp_2, h_exp_3], mourier_deep_legend{1},mourier_deep_legend{2},mourier_deep_legend{3},'Orientation','vertical');
    set(leg,'Position',[0.2 0.05 0.6 0.05]);
    leg.NumColumns = 3;
    leg.FontSize = sizeoflegendfonts;
    if flag_save_plots        
        plotname = char(compose('paper_ssr0_deep'));
        print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
        print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
        savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
        clear plotname;
    end
 clear fig; clear ax;
end
%% Paper - Mourier Comparision SSR = 0 Shallow Dimple
%==========================================================================
% Shallow dimple:
%==========================================================================
if flag.mourier_comp
    % Define figure:
    fig = figure('Units','Centimeters','Position', [0.1 0.1 34 14]);
    % Define legend:
    mourier_shallow_legend = cell(3,1);
    mourier_shallow_legend{1} = "EHL-FBNS UI";
 
    mourier_shallow_legend{2} = "EHL-FBNS QUICK";

    mourier_shallow_legend{3} = "Simulation - Mourier et al.";


    
    for plot_index = 1:length(times_to_plot)
        time_to_plot = times_to_plot(plot_index);
        subplot('Position',[0.05+(plot_index-1)*(0.15+0.05) 0.60 0.12 0.30]);
        hold on
        %==========================================================================
        % 1st Order Upwind:
        %==========================================================================
        it_sim = 2;
        % Gap height:
        plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x2;
        plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.Nx2;
        [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(ceil(plot_aux.h_ma_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

        plot_aux.h_ma_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.h_ma_line(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim));

        h1 = plot(plot_aux.h_ma_line(it_sim).x1*1e6, ...
              plot_aux.h_ma_line(it_sim).h_ma*1e9,...
              'LineWidth',widthlines,'Color',KIT_colorlist{1});
        grid on
        xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
        ylabel('$h \mathrm{[nm]}$','fontsize',sizeoffonts);  
        title(['$x_{1,d} =$ ', num2str(analysis(it_sim).reference(1+time_to_plot,3)*1e6,'%.2f'), ' $\mathrm{\mu m}$'],'fontsize',sizeoffonts)
        
        %==========================================================================
        % 2nd Order Quick:
        %==========================================================================
        it_sim = 6;
        % Gap height:
        plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x2;
        plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.Nx2;
        [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(ceil(plot_aux.h_ma_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

        plot_aux.h_ma_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.h_ma_line(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim));
        
        h2 = plot(plot_aux.h_ma_line(it_sim).x1*1e6, ...
              plot_aux.h_ma_line(it_sim).h_ma*1e9,...
              'LineWidth',widthlines,'Color',KIT_colorlist{2},'LineStyle','--');       
        %==========================================================================
        % Mourier et al. (2006):
        %==========================================================================
        % Gap height:
        plot_aux.Mourier_sim.x1 = Mourier_sim{plot_index,1}(:,1);
        plot_aux.Mourier_sim.h_ma = Mourier_sim{plot_index,1}(:,2);
        
        h3 = plot(plot_aux.Mourier_sim.x1*1e6, plot_aux.Mourier_sim.h_ma*1e9,'LineWidth', widthlines,'Color',KIT_colorlist{3},'LineStyle',':');
        set(gca, 'YDir','reverse'); 
        % Ranges: 
        xlim([-200 200]); % [micrometer]
        xticks([-200 -100 0 100 200])
        ylim([0 350]); %[nanometer]
        yticks([0 50 100 150 200 250 300 350])
        hold off
        
        % Pressures:
        subplot('Position',[0.05+(plot_index-1)*(0.15+0.05) 0.20 0.12 0.30]);
        hold on
        grid on
        %==========================================================================
        % 1st Order Upwind:
        %==========================================================================
        it_sim = 2;
        plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.hydr_pressure_line(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).sol.p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim));
        
        p1 = plot(plot_aux.hydr_pressure_line(it_sim).x1*1e6, ...
              plot_aux.hydr_pressure_line(it_sim).p_hd*1e-6,...
              'LineWidth',widthlines,'Color',KIT_colorlist{1});
        ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);   
        xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts)
        %==========================================================================
        % 2nd Order Quick:
        %==========================================================================
        it_sim = 6;
        plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.hydr_pressure_line(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).sol.p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim));
        
        p2 = plot(plot_aux.hydr_pressure_line(it_sim).x1*1e6, ...
              plot_aux.hydr_pressure_line(it_sim).p_hd*1e-6,...
              'LineWidth',widthlines,'Color',KIT_colorlist{2},'LineStyle','--');
        ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);
        %==========================================================================
        % Mourier et al. (2006):
        %==========================================================================
        plot_aux.Mourier_sim.x1 = Mourier_sim{plot_index,2}(:,1);
        plot_aux.Mourier_sim.p_hd = Mourier_sim{plot_index,2}(:,2);
        
        p3 = plot(plot_aux.Mourier_sim.x1*1e6, plot_aux.Mourier_sim.p_hd*1e-6,'LineWidth', widthlines,'Color',KIT_colorlist{3},'LineStyle',':');       
                
       % Ranges:
        xlim([-200 200]); % [micrometer]
        xticks([-200 -100 0 100 200])
        ylim([0 800]); %[MPa]
        hold off       
    end

    % Creating the legend:
    leg = legend([h1,h2,h3],mourier_shallow_legend{1},mourier_shallow_legend{2},mourier_shallow_legend{3},'Orientation','vertical');
    set(leg,'Position',[0.2 0.05 0.6 0.05]);
    leg.NumColumns = 3;
    leg.FontSize = sizeoflegendfonts;
    
    if flag_save_plots
        plotname = char(compose('paper_ssr0_shallow'));
        print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
        print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
        savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
        clear plotname;
    end  
    clear outer_row; clear p1; clear p2; clear p3; clear p_exp_1; clear p_exp_2; clear plot_index; 
    clear h1; clear h2; clear h3; clear h_exp_1; clear h_exp_2; clear h_exp_3;
end
%% Paper - Mourier Comparision SSR = -0.5 Deep Dimple
%==========================================================================
% Deep dimple:
%==========================================================================
% Set times to plot:
times_to_plot  = [177 212 231 251 305];
if flag.mourier_comp
    % Define figure:
    fig = figure('Units','Centimeters','Position', [0.1 0.1 34 14]);
    % Define legend:
    mourier_deep_legend = cell(3,1);
    mourier_deep_legend{1} = "EHL-FBNS UI";
    mourier_deep_legend{2} = "EHL-FBNS QUICK";
    mourier_deep_legend{3} = "Experiment - Mourier et al."; 

    

    for plot_index = 1:length(times_to_plot)
        time_to_plot = times_to_plot(plot_index);
        %if plot_index == 2
        %   time_to_plot = 101; %since the sampling coordinate from Mourier is not precise for the second plot a different time point is needed compared to the simulations 
        %end
        subplot('Position',[0.05+(plot_index-1)*(0.15+0.05) 0.60 0.12 0.30]);
        hold on
        % Gap height:
        %==========================================================================
        % 1st Order Upwind:
        %==========================================================================
        it_sim = 3;
        % Gap height:
        plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x2;
        plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.Nx2;
        [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(ceil(plot_aux.h_ma_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

        plot_aux.h_ma_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.h_ma_line(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim));

        h_exp_1 = plot(plot_aux.h_ma_line(it_sim).x1*1e6, ...
                   plot_aux.h_ma_line(it_sim).h_ma*1e9,...
                   'LineWidth',widthlines,'Color',KIT_colorlist{1});
        grid on
        xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
        ylabel('$h \mathrm{[nm]}$','fontsize',sizeoffonts);  
        title(['$x_{1,d} =$ ', num2str(analysis(it_sim).reference(1+time_to_plot,3)*1e6,'%.2f'), ' $\mathrm{\mu m}$'],'fontsize',sizeoffonts)
        
        %==========================================================================
        % 2nd Order Quick:
        %==========================================================================
        it_sim = 7;
        
        plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x2;
        plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.Nx2;
        [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(ceil(plot_aux.h_ma_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

        plot_aux.h_ma_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.h_ma_line(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim));
        
        h_exp_2 = plot(plot_aux.h_ma_line(it_sim).x1*1e6, ...
                   plot_aux.h_ma_line(it_sim).h_ma*1e9,...
                   'LineWidth',widthlines,'Color',KIT_colorlist{2},'LineStyle','--');
        %==========================================================================
        % Mourier et al. (2006):
        %==========================================================================
        
        plot_aux.Mourier_2_exp.x1 = Mourier_2_exp{plot_index,1}(:,1);
        plot_aux.Mourier_2_exp.h_ma = Mourier_2_exp{plot_index,1}(:,2);
        
        h_exp_3 = plot(plot_aux.Mourier_2_exp.x1*1e6, plot_aux.Mourier_2_exp.h_ma*1e9,'LineWidth', widthlines,'Color',KIT_colorlist{3},'LineStyle',':');
        
        %==========================================================================
        % Defining ranges:
        %==========================================================================
        xlim([-200 200]); % [micrometer]
        xticks([-200 -100 0 100 200])
        ylim([0 350]); %[nanometer]   
        yticks([0 50 100 150 200 250 300 350])
        set(gca, 'YDir','reverse')
        hold off
        
        % Pressures:
        subplot('Position',[0.05+(plot_index-1)*(0.15+0.05) 0.20 0.12 0.30]);
        hold on
        grid on
        %==========================================================================
        % 1st Order Upwind:
        %==========================================================================        
        it_sim = 3;
        plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.hydr_pressure_line(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).sol.p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim));
        
        p_exp_1 = plot(plot_aux.hydr_pressure_line(it_sim).x1*1e6, ...
              plot_aux.hydr_pressure_line(it_sim).p_hd*1e-6,...
              'LineWidth',widthlines,'Color',KIT_colorlist{1});
        ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts); 
        xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts)
        
        %==========================================================================
        % 2. Order QUICK
        %==========================================================================
        it_sim = 7;
        plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.hydr_pressure_line(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).sol.p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim));
        
        p_exp_2 = plot(plot_aux.hydr_pressure_line(it_sim).x1*1e6, ...
              plot_aux.hydr_pressure_line(it_sim).p_hd*1e-6,...
              'LineWidth',widthlines,'Color',KIT_colorlist{2},'LineStyle','--');
        ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);  
        % Defining ranges:
        xlim([-200 200]); % [micrometer]  
        xticks([-200 -100 0 100 200])
        ylim([0 800]); %[MPa]
        hold off
    end
    % Creating the legend:
    leg = legend([h_exp_1, h_exp_2, h_exp_3], mourier_deep_legend{1},mourier_deep_legend{2},mourier_deep_legend{3},'Orientation','vertical');
    set(leg,'Position',[0.2 0.05 0.6 0.05]);
    leg.NumColumns = 3;
    leg.FontSize = sizeoflegendfonts;
    if flag_save_plots        
        plotname = char(compose('paper_ssr-05_deep'));
        print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
        print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
        savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
        clear plotname;
    end
 clear fig; clear ax;
end
%% Paper - Mourier Comparision SSR = -0.5 Shallow Dimple
%==========================================================================
% Shallow dimple:
%==========================================================================
if flag.mourier_comp
    % Define figure:
    fig = figure('Units','Centimeters','Position', [0.1 0.1 34 14]);
    % Define legend:
    mourier_shallow_legend = cell(3,1);
    mourier_shallow_legend{1} = "EHL-FBNS UI";
    mourier_shallow_legend{2} = "EHL-FBNS QUICK";
    mourier_shallow_legend{3} = "Simulation - Mourier et al.";


    
    for plot_index = 1:length(times_to_plot)
        time_to_plot = times_to_plot(plot_index);
        subplot('Position',[0.05+(plot_index-1)*(0.15+0.05) 0.60 0.12 0.30]);
        hold on
        %==========================================================================
        % 1st Order Upwind:
        %==========================================================================
        it_sim = 4;
        % Gap height:
        plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x2;
        plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.Nx2;
        [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(ceil(plot_aux.h_ma_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

        plot_aux.h_ma_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.h_ma_line(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim));

        h1 = plot(plot_aux.h_ma_line(it_sim).x1*1e6, ...
              plot_aux.h_ma_line(it_sim).h_ma*1e9,...
              'LineWidth',widthlines,'Color',KIT_colorlist{1});
        grid on
        xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
        ylabel('$h \mathrm{[nm]}$','fontsize',sizeoffonts);  
        title(['$x_{1,d} =$ ', num2str(analysis(it_sim).reference(1+time_to_plot,3)*1e6,'%.2f'), ' $\mathrm{\mu m}$'],'fontsize',sizeoffonts)

        %==========================================================================
        % 2nd Order Quick:
        %==========================================================================
        it_sim = 8;
        % Gap height:
        plot_aux.h_ma_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x2;
        plot_aux.h_ma_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.Nx2;
        [plot_aux.h_ma_line_x2_0(it_sim),plot_aux.h_ma_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_ma_surface(it_sim).x2-plot_aux.h_ma_surface(it_sim).x2(ceil(plot_aux.h_ma_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

        plot_aux.h_ma_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.h_ma_line(it_sim).h_ma = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.h_ma(:,plot_aux.h_ma_line_x2_0_i(it_sim));
        
        h2 = plot(plot_aux.h_ma_line(it_sim).x1*1e6, ...
              plot_aux.h_ma_line(it_sim).h_ma*1e9,...
              'LineWidth',widthlines,'Color',KIT_colorlist{2},'LineStyle','--');       
        %==========================================================================
        % Mourier et al. (2006):
        %==========================================================================
        % Gap height:
        plot_aux.Mourier_2_sim.x1 = Mourier_2_sim{plot_index,1}(:,1);
        plot_aux.Mourier_2_sim.h_ma = Mourier_2_sim{plot_index,1}(:,2);
        
        h3 = plot(plot_aux.Mourier_2_sim.x1*1e6, plot_aux.Mourier_2_sim.h_ma*1e9,'LineWidth', widthlines,'Color',KIT_colorlist{3},'LineStyle',':');
        % Ranges: 
        xlim([-200 200]); % [micrometer]
        ylim([0 350]); %[nanometer]
        yticks([0 50 100 150 200 250 300 350])
        set(gca, 'YDir','reverse'); 
        hold off
        
        % Pressures:
        subplot('Position',[0.05+(plot_index-1)*(0.15+0.05) 0.20 0.12 0.30]);
        hold on
        grid on
        %==========================================================================
        % 1st Order Upwind:
        %==========================================================================
        it_sim = 4;
        plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.hydr_pressure_line(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).sol.p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim));
        
        p1 = plot(plot_aux.hydr_pressure_line(it_sim).x1*1e6, ...
              plot_aux.hydr_pressure_line(it_sim).p_hd*1e-6,...
              'LineWidth',widthlines,'Color',KIT_colorlist{1});
        ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts); 
        xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts)
        %==========================================================================
        % 2nd Order Quick:
        %==========================================================================
        it_sim = 8;
        plot_aux.hydr_pressure_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
        plot_aux.hydr_pressure_line(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).sol.p_hd(:,plot_aux.h_ma_line_x2_0_i(it_sim));
        
        p2 = plot(plot_aux.hydr_pressure_line(it_sim).x1*1e6, ...
              plot_aux.hydr_pressure_line(it_sim).p_hd*1e-6,...
              'LineWidth',widthlines,'Color',KIT_colorlist{2},'LineStyle','--');
        ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);
        %==========================================================================
        % Mourier et al. (2006):
        %==========================================================================
        plot_aux.Mourier_2_sim.x1 = Mourier_2_sim{plot_index,2}(:,1);
        plot_aux.Mourier_2_sim.p_hd = Mourier_2_sim{plot_index,2}(:,2);
        
        p3 = plot(plot_aux.Mourier_2_sim.x1*1e6, plot_aux.Mourier_2_sim.p_hd*1e-6,'LineWidth', widthlines,'Color',KIT_colorlist{3},'LineStyle',':');       
                
       % Ranges:
        xlim([-200 200]); % [micrometer]     
        xticks([-200 -100 0 100 200])
        ylim([0 800]); %[MPa]
        hold off       
    end

    % Creating the legend:
    leg = legend([h1,h2,h3],mourier_shallow_legend{1},mourier_shallow_legend{2},mourier_shallow_legend{3},'Orientation','vertical');
    set(leg,'Position',[0.2 0.05 0.6 0.05]);
    leg.NumColumns = 3;
    leg.FontSize = sizeoflegendfonts;
    if flag_save_plots
        plotname = char(compose('paper_ssr-05_shallow'));
        print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
        print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
        savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
        clear plotname;
    end  
    clear outer_row; clear p1; clear p2; clear p3; clear p_exp_1; clear p_exp_2; clear plot_index; 
    clear h1; clear h2; clear h3; clear h_exp_1; clear h_exp_2; clear h_exp_3;
end
%% Paper - Iteration Information:
fig = figure('Units','centimeters','Position',[5.1 11 figuresize(1) figuresize(2)*1.5]);
% Declare labels and colors for legend creation:
sims_to_plot_it = [1 2 3 4 5 6 7 8];
labels  = cell(length(sims_to_plot_it),1);
colors  = [1, 1, 2, 2, 3, 3, 7, 7]; % These are the indices to be picked up from the KIT colorlist
line_styles = ["-","--","-","--","-.",":","-.",":"];
plot_index = 1;

%Set axis font size:
ax = gca;
ax.FontSize = sizeofticksfonts;
hold on
for it_sim = 1:sim_N
       if ismember(it_sim,sims_to_plot_it)
           % Set legend label:
           labels{plot_index} = loaded_data(it_sim).descr;
           % Wrapper to go through all time points for this OC:
           total_time_points = loaded_data(it_sim).opc.N_t;
           for it_time = 1:total_time_points
               iteration_for_this_time    = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).alg.it_tot  ;
               x1_c_for_this_time = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).slc.x1_c;  
               plot_aux.iter(it_sim).value(it_time)         = iteration_for_this_time;
               plot_aux.iter(it_sim).dimple_center(it_time) = x1_c_for_this_time;          
           end

           plot(plot_aux.iter(it_sim).dimple_center*1e6, plot_aux.iter(it_sim).value, 'LineWidth',widthlines,'Color',KIT_colorlist{colors(plot_index)},'LineStyle',line_styles(plot_index));
           xlabel('$x_{1,d} \mathrm{[\mu m]}$','fontsize',sizeoffonts);
           ylabel('$N_n \mathrm{[-]}$','fontsize',sizeoffonts);
           xlim([-410 410])
           grid on
           ax = gca;
           ax.FontSize = sizeofticksfonts;

           plot_index = plot_index +1;
       end
end

hold off
if flag.titles_on           
   title(" Iteration information",'fontsize',sizeoffonts);   
end

lgd = legend(labels,'Location','southoutside','Orientation','horizontal');
lgd.FontSize = sizeoflegendfonts;
lgd.NumColumns = 2;

if flag_save_plots        
   plotname = char(compose('iterations_res_changed'));
   print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
   print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) ;
   print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res);
   savefig(fig,fullfile(output_main_path,append(plotname,'.fig')));
   clear plotname;
end

clear fig; clear ax; clear total_time_points; clear x1_c_for_this_time; clear iteration_for_this_time;
clear sub_result_path ;clear sub_sub_result_path; clear input_sub_sub_result_path; clear fig; clear ax;
clear labels; clear colors; clear line_styles;

% Unset Latex style:
set(groot,'defaulttextinterpreter','remove');  
set(groot,'defaultAxesTickLabelInterpreter','remove');  
set(groot,'defaultLegendInterpreter','remove');