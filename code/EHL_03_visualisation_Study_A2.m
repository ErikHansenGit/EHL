close all; clc; clearvars; 
% Visualization of the EHL solver results
%
% Recycling of script EHL_03_visualisation_Study_A for visualization of 
% grid convergence study
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
output_main_path ='./../data/EHL_03_visualisation/Study_A2/';

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
sim_N                   = 16;

% Supply simulation identifiers:
loaded_data(1).sim_id   = 'Study_A2/Nr_1/';
loaded_data(2).sim_id   = 'Study_A2/Nr_2/';
loaded_data(3).sim_id   = 'Study_A2/Nr_3/';
loaded_data(4).sim_id   = 'Study_A3/Nr_1/';
loaded_data(5).sim_id   = 'Study_A2/Nr_4/';
loaded_data(6).sim_id   = 'Study_A2/Nr_5/';
loaded_data(7).sim_id   = 'Study_A2/Nr_6/';
loaded_data(8).sim_id   = 'Study_A3/Nr_5/';
loaded_data(9).sim_id   = 'Study_A2/Nr_7/';
loaded_data(10).sim_id  = 'Study_A2/Nr_8/';
loaded_data(11).sim_id  = 'Study_A2/Nr_9/';
loaded_data(12).sim_id  = 'Study_A3/Nr_2/';
loaded_data(13).sim_id  = 'Study_A2/Nr_10/';
loaded_data(14).sim_id  = 'Study_A2/Nr_11/';
loaded_data(15).sim_id  = 'Study_A2/Nr_12/';
loaded_data(16).sim_id  = 'Study_A3/Nr_6/';

% Supply simulation descriptions:
loaded_data(1).descr    = 'Deep, UI, SSR = 0, 33';
loaded_data(2).descr    = 'Deep, UI, SSR = 0, 65';
loaded_data(3).descr    = 'Deep, UI, SSR = 0, 129';
loaded_data(4).descr    = 'Deep, UI, SSR = 0, 257';
loaded_data(5).descr    = 'Deep, QUICK, SSR = 0, 33';
loaded_data(6).descr    = 'Deep, QUICK, SSR = 0, 65';
loaded_data(7).descr    = 'Deep, QUICK, SSR = 0, 129';
loaded_data(8).descr    = 'Deep, QUICK, SSR = 0, 257';
loaded_data(9).descr    = 'Shallow, UI, SSR = 0, 33';
loaded_data(10).descr   = 'Shallow, UI, SSR = 0, 65';
loaded_data(11).descr   = 'Shallow, UI, SSR = 0, 129';
loaded_data(12).descr   = 'Shallow, UI, SSR = 0, 257';
loaded_data(13).descr   = 'Shallow, QUICK, SSR = 0, 33';
loaded_data(14).descr   = 'Shallow, QUICK, SSR = 0, 65';
loaded_data(15).descr   = 'Shallow, QUICK, SSR = 0, 129';
loaded_data(16).descr   = 'Shallow, QUICK, SSR = 0, 257';



% Choose detailed operating condition:
i_OC          = 1;

 

%==========================================================================
%% Loading the data:
%==========================================================================    
for it_sim = 1:sim_N
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


% Create a "analysis" object to reference dimple positions and times with
% time point indices:
% Initialize the object:
analysis(1).sim_id   = 'Study_A_test/Nr_1/';
analysis(2).sim_id   = 'Study_A_test/Nr_2/';
analysis(3).sim_id   = 'Study_A_test/Nr_3/';
analysis(4).sim_id   = 'Study_A_test/Nr_4/';
analysis(5).sim_id   = 'Study_A_test/Nr_5/';
analysis(6).sim_id   = 'Study_A_test/Nr_6/';

for it_sim = 1:sim_N
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

colorlist = batlowS_list(3:10);

% Size:
widthlines          = 1.7;
contour_line_width  = 1.0;
sizeoffonts         = 9;
sizeoflegendfonts   = 9;
sizeofticksfonts    = 9;
figuresize          = [13 7];


%% Convergence behaviour:
colors  = [1, 2, 3, 7]; % These are the indices to be picked up from the colorlist
linestyle = {'-','-.','--',':'};
%==========================================================================
%% Gap height:
%==========================================================================
% Deep dimple:
%==========================================================================
% 1st Order Upwind:
%==========================================================================

% Define figure:
fig = figure('Units','Centimeters','Position', [2 20 figuresize]);
for plot_index = 1:4
    it_sim = plot_index;

    time_to_plot = ceil(loaded_data(it_sim).opc.N_t/2);

    % Gap height:
    plot_aux.h_m_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x2;
    plot_aux.h_m_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.Nx2;
    [plot_aux.h_m_line_x2_0(it_sim),plot_aux.h_m_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_m_surface(it_sim).x2-plot_aux.h_m_surface(it_sim).x2(ceil(plot_aux.h_m_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

    plot_aux.h_m_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
    plot_aux.h_m_line(it_sim).h_m = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.h_m(:,plot_aux.h_m_line_x2_0_i(it_sim));

    plot(plot_aux.h_m_line(it_sim).x1*1e6, ...
               plot_aux.h_m_line(it_sim).h_m*1e9,...
               'LineWidth',widthlines,'Color',colorlist{colors(plot_index)},'linestyle',linestyle{plot_index});
    grid on
    xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    ylabel('$h \mathrm{[nm]}$','fontsize',sizeoffonts);  
    set(gca, 'YDir','reverse'); 
    %==========================================================================
    % Defining ranges:
    %==========================================================================
    xlim([-200 200]); % [micrometer]
    ylim([0 350]); %[nanometer]        
    hold on
end
legend({'$33^3$','$65^3$','$129^3$','$257^3$'},'fontsize',sizeoffonts)
hold off


if flag_save_plots        
    plotname = char(compose('paper_ssr0_deep_UI_h_convergence'));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;

%==========================================================================
% 2nd Order Upwind:
%==========================================================================

% Define figure:
fig = figure('Units','Centimeters','Position', [20 20 figuresize]);
for plot_index = 5:8
    it_sim = plot_index;

    time_to_plot = ceil(loaded_data(it_sim).opc.N_t/2);

    % Gap height:
    plot_aux.h_m_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x2;
    plot_aux.h_m_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.Nx2;
    [plot_aux.h_m_line_x2_0(it_sim),plot_aux.h_m_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_m_surface(it_sim).x2-plot_aux.h_m_surface(it_sim).x2(ceil(plot_aux.h_m_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

    plot_aux.h_m_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
    plot_aux.h_m_line(it_sim).h_m = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.h_m(:,plot_aux.h_m_line_x2_0_i(it_sim));

    plot(plot_aux.h_m_line(it_sim).x1*1e6, ...
               plot_aux.h_m_line(it_sim).h_m*1e9,...
               'LineWidth',widthlines,'Color',colorlist{colors(plot_index-4)},'linestyle',linestyle{plot_index-4});
    grid on
    xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    ylabel('$h \mathrm{[nm]}$','fontsize',sizeoffonts);  
    set(gca, 'YDir','reverse'); 
    %==========================================================================
    % Defining ranges:
    %==========================================================================
    xlim([-200 200]); % [micrometer]
    ylim([0 350]); %[nanometer]        
    hold on
end
legend({'$33^3$','$65^3$','$129^3$','$257^3$'},'fontsize',sizeoffonts)
hold off


if flag_save_plots        
    plotname = char(compose('paper_ssr0_deep_QUICK_h_convergence'));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;

%==========================================================================
% Shallow dimple:
%==========================================================================
% 1st Order Upwind:
%==========================================================================

% Define figure:
fig = figure('Units','Centimeters','Position', [6 20 figuresize]);
for plot_index = 9:12
    it_sim = plot_index;

    time_to_plot = ceil(loaded_data(it_sim).opc.N_t/2);

    % Gap height:
    plot_aux.h_m_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x2;
    plot_aux.h_m_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.Nx2;
    [plot_aux.h_m_line_x2_0(it_sim),plot_aux.h_m_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_m_surface(it_sim).x2-plot_aux.h_m_surface(it_sim).x2(ceil(plot_aux.h_m_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

    plot_aux.h_m_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
    plot_aux.h_m_line(it_sim).h_m = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.h_m(:,plot_aux.h_m_line_x2_0_i(it_sim));

    plot(plot_aux.h_m_line(it_sim).x1*1e6, ...
               plot_aux.h_m_line(it_sim).h_m*1e9,...
               'LineWidth',widthlines,'Color',colorlist{colors(plot_index-8)},'linestyle',linestyle{plot_index-8});
    grid on
    xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    ylabel('$h \mathrm{[nm]}$','fontsize',sizeoffonts);  
    set(gca, 'YDir','reverse'); 
    %==========================================================================
    % Defining ranges:
    %==========================================================================
    xlim([-200 200]); % [micrometer]
    ylim([0 350]); %[nanometer]        
    hold on
end
legend({'$33^3$','$65^3$','$129^3$','$257^3$'},'fontsize',sizeoffonts)
hold off


if flag_save_plots        
    plotname = char(compose('paper_ssr0_shallow_UI_h_convergence'));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;

%==========================================================================
% 2nd Order Upwind:
%==========================================================================

% Define figure:
fig = figure('Units','Centimeters','Position', [24 20 figuresize]);
for plot_index = 13:16
    it_sim = plot_index;

    time_to_plot = ceil(loaded_data(it_sim).opc.N_t/2);

    % Gap height:
    plot_aux.h_m_surface(it_sim).x2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x2;
    plot_aux.h_m_surface(it_sim).Nx2 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.Nx2;
    [plot_aux.h_m_line_x2_0(it_sim),plot_aux.h_m_line_x2_0_i(it_sim)] = min(abs(plot_aux.h_m_surface(it_sim).x2-plot_aux.h_m_surface(it_sim).x2(ceil(plot_aux.h_m_surface(it_sim).Nx2/2))));       % [m]   x2-position of line plot

    plot_aux.h_m_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
    plot_aux.h_m_line(it_sim).h_m = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.h_m(:,plot_aux.h_m_line_x2_0_i(it_sim));

    plot(plot_aux.h_m_line(it_sim).x1*1e6, ...
               plot_aux.h_m_line(it_sim).h_m*1e9,...
               'LineWidth',widthlines,'Color',colorlist{colors(plot_index-12)},'linestyle',linestyle{plot_index-12});
    grid on
    xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    ylabel('$h \mathrm{[nm]}$','fontsize',sizeoffonts);  
    set(gca, 'YDir','reverse'); 
    %==========================================================================
    % Defining ranges:
    %==========================================================================
    xlim([-200 200]); % [micrometer]
    ylim([0 350]); %[nanometer]        
    hold on
end
legend({'$33^3$','$65^3$','$129^3$','$257^3$'},'fontsize',sizeoffonts)
hold off


if flag_save_plots        
    plotname = char(compose('paper_ssr0_shallow_QUICK_h_convergence'));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;

%==========================================================================
%% Central height:
%==========================================================================
% Deep dimple:
%==========================================================================
% 1st Order Upwind:
%==========================================================================
fig = figure('Units','Centimeters','Position', [2 1 figuresize]);
for plot_index = 1:4
    it_sim = plot_index;
    plot_aux.h_C_line(it_sim).x1_c = zeros(loaded_data(it_sim).opc.N_t,1);
    plot_aux.h_C_line(it_sim).h_C = zeros(loaded_data(it_sim).opc.N_t,1);
    
    plot_aux.h_m_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(1).h.x1;
    plot_aux.h_m_surface(it_sim).Nx1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(1).h.Nx1;
    [plot_aux.h_m_line_x1_0(it_sim),plot_aux.h_m_line_x1_0_i(it_sim)] = min(abs(plot_aux.h_m_surface(it_sim).x1-plot_aux.h_m_surface(it_sim).x1(ceil(plot_aux.h_m_surface(it_sim).Nx1/2))));       % [m]   x2-position of line plot

    
    for it_time = 1:loaded_data(it_sim).opc.N_t
        plot_aux.h_C_line(it_sim).x1_c(it_time) = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).slc.x1_c;
        plot_aux.h_C_line(it_sim).h_C(it_time) = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.h_m(plot_aux.h_m_line_x1_0_i(it_sim),plot_aux.h_m_line_x2_0_i(it_sim));
    end

    plot(plot_aux.h_C_line(it_sim).x1_c*1e6, ...
               plot_aux.h_C_line(it_sim).h_C*1e9,...
               'LineWidth',widthlines,'Color',colorlist{colors(plot_index)},'linestyle',linestyle{plot_index});
    grid on
    xlabel('$x_{1,d} \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    ylabel('$h_C \mathrm{[nm]}$','fontsize',sizeoffonts);  
    set(gca, 'YDir','reverse'); 
    %==========================================================================
    % Defining ranges:
    %==========================================================================
    xlim([-200 200]); % [micrometer]
    ylim([0 350]); %[nanometer]        
    hold on
end
legend({'$33^3$','$65^3$','$129^3$','$257^3$'},'fontsize',sizeoffonts)
hold off


if flag_save_plots        
    plotname = char(compose('paper_ssr0_deep_UI_h_C_convergence'));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;

%==========================================================================
% 2nd Order Upwind:
%==========================================================================

fig = figure('Units','Centimeters','Position', [20 1 figuresize]);
for plot_index = 5:8
    it_sim = plot_index;
    plot_aux.h_C_line(it_sim).x1_c = zeros(loaded_data(it_sim).opc.N_t,1);
    plot_aux.h_C_line(it_sim).h_C = zeros(loaded_data(it_sim).opc.N_t,1);
    
    plot_aux.h_m_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(1).h.x1;
    plot_aux.h_m_surface(it_sim).Nx1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(1).h.Nx1;
    [plot_aux.h_m_line_x1_0(it_sim),plot_aux.h_m_line_x1_0_i(it_sim)] = min(abs(plot_aux.h_m_surface(it_sim).x1-plot_aux.h_m_surface(it_sim).x1(ceil(plot_aux.h_m_surface(it_sim).Nx1/2))));       % [m]   x2-position of line plot

    
    for it_time = 1:loaded_data(it_sim).opc.N_t
        plot_aux.h_C_line(it_sim).x1_c(it_time) = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).slc.x1_c;
        plot_aux.h_C_line(it_sim).h_C(it_time) = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.h_m(plot_aux.h_m_line_x1_0_i(it_sim),plot_aux.h_m_line_x2_0_i(it_sim));
    end

    plot(plot_aux.h_C_line(it_sim).x1_c*1e6, ...
               plot_aux.h_C_line(it_sim).h_C*1e9,...
               'LineWidth',widthlines,'Color',colorlist{colors(plot_index-4)},'linestyle',linestyle{plot_index-4});
    grid on
    xlabel('$x_{1,d} \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    ylabel('$h_C \mathrm{[nm]}$','fontsize',sizeoffonts);  
    set(gca, 'YDir','reverse'); 
    %==========================================================================
    % Defining ranges:
    %==========================================================================
    xlim([-200 200]); % [micrometer]
    ylim([0 350]); %[nanometer]        
    hold on
end
legend({'$33^3$','$65^3$','$129^3$','$257^3$'},'fontsize',sizeoffonts)
hold off


if flag_save_plots        
    plotname = char(compose('paper_ssr0_deep_QUICK_h_C_convergence'));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;

%==========================================================================
% Shallow dimple:
%==========================================================================
% 1st Order Upwind:
%==========================================================================
fig = figure('Units','Centimeters','Position', [6 1 figuresize]);
for plot_index = 9:12
    it_sim = plot_index;
    plot_aux.h_C_line(it_sim).x1_c = zeros(loaded_data(it_sim).opc.N_t,1);
    plot_aux.h_C_line(it_sim).h_C = zeros(loaded_data(it_sim).opc.N_t,1);
    
    plot_aux.h_m_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(1).h.x1;
    plot_aux.h_m_surface(it_sim).Nx1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(1).h.Nx1;
    [plot_aux.h_m_line_x1_0(it_sim),plot_aux.h_m_line_x1_0_i(it_sim)] = min(abs(plot_aux.h_m_surface(it_sim).x1-plot_aux.h_m_surface(it_sim).x1(ceil(plot_aux.h_m_surface(it_sim).Nx1/2))));       % [m]   x2-position of line plot

    
    for it_time = 1:loaded_data(it_sim).opc.N_t
        plot_aux.h_C_line(it_sim).x1_c(it_time) = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).slc.x1_c;
        plot_aux.h_C_line(it_sim).h_C(it_time) = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.h_m(plot_aux.h_m_line_x1_0_i(it_sim),plot_aux.h_m_line_x2_0_i(it_sim));
    end

    plot(plot_aux.h_C_line(it_sim).x1_c*1e6, ...
               plot_aux.h_C_line(it_sim).h_C*1e9,...
               'LineWidth',widthlines,'Color',colorlist{colors(plot_index-8)},'linestyle',linestyle{plot_index-8});
    grid on
    xlabel('$x_{1,d} \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    ylabel('$h_C \mathrm{[nm]}$','fontsize',sizeoffonts);  
    set(gca, 'YDir','reverse'); 
    %==========================================================================
    % Defining ranges:
    %==========================================================================
    xlim([-200 200]); % [micrometer]
    ylim([0 350]); %[nanometer]        
    hold on
end
legend({'$33^3$','$65^3$','$129^3$','$257^3$'},'fontsize',sizeoffonts)
hold off


if flag_save_plots        
    plotname = char(compose('paper_ssr0_shallow_UI_h_C_convergence'));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;

%==========================================================================
% 2nd Order Upwind:
%==========================================================================

fig = figure('Units','Centimeters','Position', [24 1 figuresize]);
for plot_index = 13:16
    it_sim = plot_index;
    plot_aux.h_C_line(it_sim).x1_c = zeros(loaded_data(it_sim).opc.N_t,1);
    plot_aux.h_C_line(it_sim).h_C = zeros(loaded_data(it_sim).opc.N_t,1);
    
    plot_aux.h_m_surface(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(1).h.x1;
    plot_aux.h_m_surface(it_sim).Nx1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(1).h.Nx1;
    [plot_aux.h_m_line_x1_0(it_sim),plot_aux.h_m_line_x1_0_i(it_sim)] = min(abs(plot_aux.h_m_surface(it_sim).x1-plot_aux.h_m_surface(it_sim).x1(ceil(plot_aux.h_m_surface(it_sim).Nx1/2))));       % [m]   x2-position of line plot

    
    for it_time = 1:loaded_data(it_sim).opc.N_t
        plot_aux.h_C_line(it_sim).x1_c(it_time) = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).slc.x1_c;
        plot_aux.h_C_line(it_sim).h_C(it_time) = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(it_time).h.h_m(plot_aux.h_m_line_x1_0_i(it_sim),plot_aux.h_m_line_x2_0_i(it_sim));
    end

    plot(plot_aux.h_C_line(it_sim).x1_c*1e6, ...
               plot_aux.h_C_line(it_sim).h_C*1e9,...
               'LineWidth',widthlines,'Color',colorlist{colors(plot_index-12)},'linestyle',linestyle{plot_index-12});
    grid on
    xlabel('$x_{1,d} \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    ylabel('$h_C \mathrm{[nm]}$','fontsize',sizeoffonts);  
    set(gca, 'YDir','reverse'); 
    %==========================================================================
    % Defining ranges:
    %==========================================================================
    xlim([-200 200]); % [micrometer]
    ylim([0 350]); %[nanometer]        
    hold on
end
legend({'$33^3$','$65^3$','$129^3$','$257^3$'},'fontsize',sizeoffonts)
hold off


if flag_save_plots        
    plotname = char(compose('paper_ssr0_shallow_QUICK_h_C_convergence'));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;

%==========================================================================
%% Pressure:
%==========================================================================
% Deep dimple:
%==========================================================================
% 1st Order Upwind:
%==========================================================================

% Define figure:
fig = figure('Units','Centimeters','Position', [2 10 figuresize]);
for plot_index = 1:4
    it_sim = plot_index;

    time_to_plot = ceil(loaded_data(it_sim).opc.N_t/2);

    % Pressure:
    plot_aux.p_hd_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
    plot_aux.p_hd_line(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).sol.p_hd(:,plot_aux.h_m_line_x2_0_i(it_sim));

    plot(plot_aux.h_m_line(it_sim).x1*1e6, ...
               plot_aux.p_hd_line(it_sim).p_hd*1e-6,...
               'LineWidth',widthlines,'Color',colorlist{colors(plot_index)},'linestyle',linestyle{plot_index});
    grid on
    xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);  
    %==========================================================================
    % Defining ranges:
    %==========================================================================
    xlim([-200 200]); % [micrometer]
    ylim([0 1000]); % [MPa]        
    hold on
end
legend({'$33^3$','$65^3$','$129^3$','$257^3$'},'fontsize',sizeoffonts)
hold off


if flag_save_plots        
    plotname = char(compose('paper_ssr0_deep_UI_p_hd_convergence'));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;

%==========================================================================
% 2nd Order Upwind:
%==========================================================================

% Define figure:
fig = figure('Units','Centimeters','Position', [20 10 figuresize]);
for plot_index = 5:8
    it_sim = plot_index;

    time_to_plot = ceil(loaded_data(it_sim).opc.N_t/2);

    % Pressure:
    plot_aux.p_hd_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
    plot_aux.p_hd_line(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).sol.p_hd(:,plot_aux.h_m_line_x2_0_i(it_sim));

    plot(plot_aux.h_m_line(it_sim).x1*1e6, ...
               plot_aux.p_hd_line(it_sim).p_hd*1e-6,...
               'LineWidth',widthlines,'Color',colorlist{colors(plot_index-4)},'linestyle',linestyle{plot_index-4});
    grid on
    xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);  
    %==========================================================================
    % Defining ranges:
    %==========================================================================
    xlim([-200 200]); % [micrometer]
    ylim([0 1000]); % [MPa]        
    hold on
end
legend({'$33^3$','$65^3$','$129^3$','$257^3$'},'fontsize',sizeoffonts)
hold off


if flag_save_plots        
    plotname = char(compose('paper_ssr0_deep_QUICK_p_hd_convergence'));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;

%==========================================================================
% Shallow dimple:
%==========================================================================
% 1st Order Upwind:
%==========================================================================

% Define figure:
fig = figure('Units','Centimeters','Position', [6 10 figuresize]);
for plot_index = 9:12
    it_sim = plot_index;

    time_to_plot = ceil(loaded_data(it_sim).opc.N_t/2);

    % Pressure:
    plot_aux.p_hd_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
    plot_aux.p_hd_line(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).sol.p_hd(:,plot_aux.h_m_line_x2_0_i(it_sim));

    plot(plot_aux.h_m_line(it_sim).x1*1e6, ...
               plot_aux.p_hd_line(it_sim).p_hd*1e-6,...
               'LineWidth',widthlines,'Color',colorlist{colors(plot_index-8)},'linestyle',linestyle{plot_index-8});
    grid on
    xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);  
    %==========================================================================
    % Defining ranges:
    %==========================================================================
    xlim([-200 200]); % [micrometer]
    ylim([0 1000]); % [MPa]        
    hold on
end
legend({'$33^3$','$65^3$','$129^3$','$257^3$'},'fontsize',sizeoffonts)
hold off


if flag_save_plots        
    plotname = char(compose('paper_ssr0_shallow_UI_p_hd_convergence'));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;

%==========================================================================
% 2nd Order Upwind:
%==========================================================================

% Define figure:
fig = figure('Units','Centimeters','Position', [24 10 figuresize]);
for plot_index = 13:16
    it_sim = plot_index;

    time_to_plot = ceil(loaded_data(it_sim).opc.N_t/2);

    % Pressure:
    plot_aux.p_hd_line(it_sim).x1 = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).h.x1;
    plot_aux.p_hd_line(it_sim).p_hd = loaded_data(it_sim).loaded_OC(i_OC).loaded_Time(time_to_plot).sol.p_hd(:,plot_aux.h_m_line_x2_0_i(it_sim));

    plot(plot_aux.h_m_line(it_sim).x1*1e6, ...
               plot_aux.p_hd_line(it_sim).p_hd*1e-6,...
               'LineWidth',widthlines,'Color',colorlist{colors(plot_index-12)},'linestyle',linestyle{plot_index-12});
    grid on
    xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
    ylabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts);  
    %==========================================================================
    % Defining ranges:
    %==========================================================================
    xlim([-200 200]); % [micrometer]
    ylim([0 1000]); % [MPa]        
    hold on
end
legend({'$33^3$','$65^3$','$129^3$','$257^3$'},'fontsize',sizeoffonts)
hold off


if flag_save_plots        
    plotname = char(compose('paper_ssr0_shallow_QUICK_p_hd_convergence'));
    print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
    print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) 
    print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res)
    savefig(fig,fullfile(output_main_path,append(plotname,'.fig')))
    clear plotname;
end
clear fig; clear ax;

%% Iteration Information:
% Deep dimple
fig = figure('Units','centimeters','Position',[5.1 11 figuresize(1) figuresize(2)*1.5]);
% Declare labels and colors for legend creation:
% sims_to_plot_it = [1 2 3 4 5 6 7 8];
% labels  = cell(length(sims_to_plot_it),1);
colors  = [1, 2, 3, 7, 1, 2, 3, 7]; % These are the indices to be picked up from the KIT colorlist
line_styles = ["-","-","-","-","--","--","--","--"];
plot_index = 1;

%Set axis font size:
ax = gca;
ax.FontSize = sizeofticksfonts;
hold on
for it_sim = 1:8
           % Set legend label:
%            labels{plot_index} = loaded_data(it_sim).descr;
           % Wrapper to go through all time points for this OC:
           total_time_points = loaded_data(it_sim).opc.N_t;
           for it_time = 1:total_time_points
               iteration_for_this_time    = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).alg.it_tot  ;
               x1_c_for_this_time = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).slc.x1_c;  
               plot_aux.iter(it_sim).value(it_time)         = iteration_for_this_time;
               plot_aux.iter(it_sim).dimple_center(it_time) = x1_c_for_this_time;          
           end
           %%%mod!!
           plot(plot_aux.iter(it_sim).dimple_center*1e6, plot_aux.iter(it_sim).value, 'LineWidth',widthlines,'Color',colorlist{colors(plot_index)},'LineStyle',line_styles(plot_index));
           xlabel('$x_{1,d} \mathrm{[\mu m]}$','fontsize',sizeoffonts);
           ylabel('$N_n \mathrm{[-]}$','fontsize',sizeoffonts);
           xlim([-410 410])
           grid on
           ax = gca;
           ax.FontSize = sizeofticksfonts;

           plot_index = plot_index +1;
end

hold off
% if flag.titles_on           
%    title(" Iteration information",'fontsize',sizeoffonts);   
% end
lgd = legend({'UI $33^3$','UI $65^3$','UI $129^3$','UI $257^3$','QUICK $33^3$','QUICK $65^3$','QUICK $129^3$','QUICK $257^3$'},'Location','southoutside','Orientation','horizontal','fontsize',sizeoffonts);
lgd.FontSize = sizeoflegendfonts;
lgd.NumColumns = 2;

if flag_save_plots        
   plotname = char(compose('iterations_deep'));
   print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
   print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) ;
   print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res);
   savefig(fig,fullfile(output_main_path,append(plotname,'.fig')));
   clear plotname;
end

clear fig; clear ax; 

% Shallow dimple
fig = figure('Units','centimeters','Position',[5.1 11 figuresize(1) figuresize(2)*1.5]);
% Declare labels and colors for legend creation:
% sims_to_plot_it = 9:16;
% labels  = cell(length(sims_to_plot_it),1);
colors  = [1, 2, 3, 7, 1, 2, 3, 7]; % These are the indices to be picked up from the KIT colorlist
line_styles = ["-","-","-","-","--","--","--","--"];

%Set axis font size:
ax = gca;
ax.FontSize = sizeofticksfonts;
hold on
for it_sim = 9:16
           % Set legend label:
%            labels{plot_index} = loaded_data(it_sim).descr;
           % Wrapper to go through all time points for this OC:
           total_time_points = loaded_data(it_sim).opc.N_t;
           for it_time = 1:total_time_points
               iteration_for_this_time    = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).alg.it_tot  ;
               x1_c_for_this_time = loaded_data(it_sim).loaded_OC.loaded_Time(it_time).slc.x1_c;  
               plot_aux.iter(it_sim).value(it_time)         = iteration_for_this_time;
               plot_aux.iter(it_sim).dimple_center(it_time) = x1_c_for_this_time;          
           end

           plot(plot_aux.iter(it_sim).dimple_center*1e6, plot_aux.iter(it_sim).value, 'LineWidth',widthlines,'Color',colorlist{colors(plot_index-8)},'LineStyle',line_styles(plot_index-8));
           xlabel('$x_{1,d} \mathrm{[\mu m]}$','fontsize',sizeoffonts);
           ylabel('$N_n \mathrm{[-]}$','fontsize',sizeoffonts);
           xlim([-410 410])
           grid on
           ax = gca;
           ax.FontSize = sizeofticksfonts;

           plot_index = plot_index +1;
end

hold off
% if flag.titles_on           
%    title(" Iteration information",'fontsize',sizeoffonts);   
% end

lgd = legend({'UI $33^3$','UI $65^3$','UI $129^3$','UI $257^3$','QUICK $33^3$','QUICK $65^3$','QUICK $129^3$','QUICK $257^3$'},'Location','southoutside','Orientation','horizontal','fontsize',sizeoffonts);
lgd.FontSize = sizeoflegendfonts;
lgd.NumColumns = 2;

if flag_save_plots        
   plotname = char(compose('iterations_shallow'));
   print(fig,fullfile(output_main_path,append(plotname,'.eps')),'-depsc',print_res)
   print(fig,fullfile(output_main_path,append(plotname,'.svg')),'-dsvg',print_res) ;
   print(fig,fullfile(output_main_path,append(plotname,'.png')),'-dpng',print_res);
   savefig(fig,fullfile(output_main_path,append(plotname,'.fig')));
   clear plotname;
end

clear fig; clear ax; clear total_time_points; clear x1_c_for_this_time; clear iteration_for_this_time;
clear sub_result_path ;clear sub_sub_result_path; clear input_sub_sub_result_path; clear fig; clear ax;
clear labels; clear colors; clear line_styles;

%==========================================================================
% Unset Latex style:
set(groot,'defaulttextinterpreter','remove');  
set(groot,'defaultAxesTickLabelInterpreter','remove');  
set(groot,'defaultLegendInterpreter','remove');