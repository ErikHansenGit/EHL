close all; clc; clear all; 
% Visualization of the EHL solver results
%
% This script imports the data exported by the script
% "EHL_02_mainprocess.m" and visualizes some of its contents in Latex style.
% The input properties used in the solver "EHL_02_mainprocess.m" are 
% consolidated in a table which is specified in the "Write table" section of this script.
% The exact file paths of the in- and output data can be
% specified in the "File path information" section of this script. The
% evaluated operating condition, the boolean whether to save plots and
% table and the printing resolution are also specified in this section.
% The plot settings can be adjusted in the "Plot settings" section of this
% script.
% 
% Erik Hansen, 07.09.2020

%% File path information:
% Input path:
input_main_path = sprintf('%s','./../data/Rosenkranz/EHL_02_mainprocess/Output');

% Output path:
output_main_path = sprintf('%s','./../data/Rosenkranz/EHL_03_visualization/Output');
% Save results:
output_result_path = fullfile(output_main_path,'/Plots');
output_result_table_path = fullfile(output_main_path,'/Tables');

% Choose operating conditon number for detailed plot:
i_OC                = 1;                                        % [-]

% Choose operating conditon number for detailed plot:
flag_save_plots     = false;                                    % [-]   boolean whether to save the plots and table or not
print_res           = '-r600';                                  % [-]   resolution used when the plots are printed

% Create output directiories:
if flag_save_plots 
    mkdir (output_main_path)
    mkdir (output_result_path)
    mkdir (output_result_table_path)
end

%% Load input information:
input_used_input_path = fullfile(input_main_path,'/Used_input');
load(fullfile(input_used_input_path,'/fld.mat'));
load(fullfile(input_used_input_path,'/sld.mat'));
load(fullfile(input_used_input_path,'/opc.mat'));
clear input_used_input_path;
% Load result information:
input_result_path = fullfile(input_main_path,'/Result');
load(fullfile(input_result_path,'/str.mat'));
% Load detailed result information:
sub_result_path = sprintf('/OC_%i',i_OC);
input_sub_result_path = fullfile(input_result_path,sub_result_path);
clear input_result_path; clear sub_result_path;
load(fullfile(input_sub_result_path,'/alg.mat'));
load(fullfile(input_sub_result_path,'/h.mat'));
load(fullfile(input_sub_result_path,'/prop.mat'));
load(fullfile(input_sub_result_path,'/sol.mat'));
load(fullfile(input_sub_result_path,'/res.mat'));
clear input_sub_result_path;

%% Plot settings:
% Line plots:
KIT_colorlist={[0,150,130]/255,[162 34 35]/255,[70 100 170]/255,[252 229 0]/255,[140 182 60]/256,[223 155 27]/255,[167 130 46]/255,[163 16 124]/255,[35 161 224]/255};
% Surface plots:
colmap_h = parula;
colmap_p = parula;

% Size:
widthlines          = 1.2;
contour_line_width  = 1.2;
sizeoffonts         = 11;
sizeoflegendfonts   = 11;
sizeofticksfonts    = 7;
figuresize          = [7 7];
% Ranges and ticks:
U_range             = [1e-1 1];                                             % [m/s]
U_ticks             = opc.u_up - opc.u_low;                                 % [m/s]
C_f_range           = [1e-2 1e0];                                           % [-]
C_f_ticks           = logspace(log10(C_f_range(1)),log10(C_f_range(2)),3);  % [-]
h_range             = [0 500];                                              % [nm]
h_ticks             = linspace(h_range(1),h_range(2),11);                   % [nm]
h_cont_range        = [0 500];                                              % [nm]
h_cont_ticks        = linspace(h_cont_range(1),h_cont_range(2),50);         % [nm]
p_range             = [0 1500];                                             % [MPa]
p_ticks             = linspace(p_range(1),p_range(2),11);                   % [MPa]
tau_range           = [0 200];                                              % [MPa]
tau_ticks           = linspace(tau_range(1),tau_range(2),11);               % [MPa]
t_range             = [0 1];                                                % [-]
t_ticks             = linspace(t_range(1),t_range(2),6);                    % [-]
x_n_ticks           = 11;                                                   % [-]
x1_range            = [-100 100];                                           % [\mu m]
x1_ticks            = linspace(x1_range(1),x1_range(2),x_n_ticks);          % [\mu m]
x2_range            = [-100 100];                                           % [\mu m]
x2_ticks            = linspace(x2_range(1),x2_range(2),x_n_ticks);          % [\mu m]

%% Plot:
% Latex style:
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

% General plots:
% Stribeck curve:
fig = figure('Units','centimeters','Position',[6.1 11 figuresize]);
loglog(opc.u_up-opc.u_low,str.C_f,'LineWidth',widthlines,'Color',KIT_colorlist{1})
grid on
ax = gca;
ax.FontSize = sizeofticksfonts;
xlim(U_range)
xticks(U_ticks)
ylim(C_f_range)
yticks(C_f_ticks)
xlabel('$U_2-U_1 \mathrm{[m/s]}$','fontsize',sizeoffonts);
ylabel('$C_f \mathrm{[-]}$','fontsize',sizeoffonts);
if flag_save_plots 
    print(fig,fullfile(output_result_path,'/C_f.svg'),'-dsvg',print_res) 
    print(fig,fullfile(output_result_path,'/C_f.png'),'-dpng',print_res)
    savefig(fig,fullfile(output_result_path,'/C_f.fig'))
end
clear fig; clear ax;

% Minimum gap height:
fig = figure('Units','centimeters','Position',[24.1 11 figuresize]);
loglog(opc.u_up-opc.u_low,str.h_ma_ad_min*1e9,'LineWidth',widthlines,'Color',KIT_colorlist{1})
grid on
ax = gca;
ax.FontSize = sizeofticksfonts;
xlim(U_range)
xticks(U_ticks)
ylim(h_range)
yticks(h_ticks)
xlabel('$U_2-U_1 \mathrm{[m/s]}$','fontsize',sizeoffonts);
ylabel('$min(h_{ma,ad}) \mathrm{[nm]}$','fontsize',sizeoffonts);
if flag_save_plots 
    print(fig,fullfile(output_result_path,'/h_ad_min.svg'),'-dsvg',print_res) 
    print(fig,fullfile(output_result_path,'/h_ad_min.png'),'-dpng',print_res)
    savefig(fig,fullfile(output_result_path,'/h_ad_min.fig'))
end
clear fig; clear ax;

% Number of iterations:
fig = figure('Units','centimeters','Position',[18.1 11 figuresize]);
loglog(opc.u_up-opc.u_low,str.it_tot,'LineWidth',widthlines,'Color',KIT_colorlist{1})
hold on
loglog(opc.u_up,str.it_W,'LineWidth',widthlines,'Color',KIT_colorlist{2})
grid on
ax = gca;
ax.FontSize = sizeofticksfonts;
xlim(U_range)
xticks(U_ticks)
xlabel('$U_2-U_1 \mathrm{[m/s]}$','fontsize',sizeoffonts);
ylabel('$N_{it} \mathrm{[-]}$','fontsize',sizeoffonts);
lgd = legend('$N_{it,tot}$','$N_{it,W}$');
lgd.FontSize = sizeoflegendfonts;
if flag_save_plots 
    print(fig,fullfile(output_result_path,'/N_it.svg'),'-dsvg',print_res) 
    print(fig,fullfile(output_result_path,'/N_it.png'),'-dpng',print_res)
    savefig(fig,fullfile(output_result_path,'/N_it.fig'))
end
clear fig; clear ldg; clear ax;

% Detailed plots:
fprintf('\n-----------------------------\n');
fprintf('\nOC_i     = %d\n',i_OC);
fprintf('\nU        = %d m/s\n',opc.u_up(i_OC));
fprintf('\n-----------------------------\n');

% Residuals:
fig = figure('Units','centimeters','Position',[12.1 11 figuresize]);
it_tot = 1:alg.it_tot;
semilogy(it_tot,res.delta_p,'-','LineWidth',widthlines,'Color',KIT_colorlist{1})
hold on
semilogy(it_tot,res.delta_t,'-.','LineWidth',widthlines,'Color',KIT_colorlist{2})
hold on
semilogy(res.W(:,1),res.W(:,2),'*','LineWidth',widthlines,'Color',KIT_colorlist{3})
grid on
ax = gca;
ax.FontSize = sizeofticksfonts;
xlabel('$it_{tot} \mathrm{[-]}$','fontsize',sizeoffonts);
ylabel('$res \mathrm{[-]}$','fontsize',sizeoffonts);
lgd = legend('$r_{\delta p,hd}$','$r_{\theta}$','$r_{W}$');
lgd.FontSize = sizeoflegendfonts;
if flag_save_plots 
    print(fig,fullfile(output_result_path,'/residuals.svg'),'-dsvg',print_res) 
    print(fig,fullfile(output_result_path,'/residuals.png'),'-dpng',print_res)
    savefig(fig,fullfile(output_result_path,'/residuals.fig'))
end
clear fig; clear lgd;
clear it_tot; clear ax;

% Contact pressure:
fig = figure('Units','centimeters','Position',[6.1 01 figuresize]);
surf(h.x1*1e6,h.x2*1e6,sol.p_con'*1e-6)
shading interp
material dull
colormap(colmap_p)
camlight
ax = gca;
ax.FontSize = sizeofticksfonts;
xlim(x1_range)
xticks(x1_ticks)
ylim(x2_range)
yticks(x2_ticks)
zlim(p_range)
zticks(p_ticks)
xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
zlabel('$p_{con} \mathrm{[MPa]}$','fontsize',sizeoffonts)
if flag_save_plots 
    print(fig,fullfile(output_result_path,'/p_con.svg'),'-dsvg',print_res) 
    print(fig,fullfile(output_result_path,'/p_con.png'),'-dpng',print_res)
    savefig(fig,fullfile(output_result_path,'/p_con.fig'))
end
clear fig; clear ax;

% Fluid pressure:
% Line along x1:
[x2_0,x2_0_i] = min(abs(h.x2));
fig = figure('Units','centimeters','Position',[12.1 01 figuresize]);
plot(h.x1*1e6,sol.p_hd(:,x2_0_i)*1e-6,'LineWidth',widthlines,'Color',KIT_colorlist{1})
grid on
ax = gca;
ax.FontSize = sizeofticksfonts;
xlim(x1_range)
xticks(x1_ticks)
ylim(p_range)
yticks(p_ticks)
xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
ylabel('$p_hd \mathrm{[MPa]}$','fontsize',sizeoffonts);
title("$x_2=$ " + x2_0*1e6 + "$\mathrm{\mu m}$",'fontsize',sizeoffonts)
if flag_save_plots 
    print(fig,fullfile(output_result_path,'/p_hd_x1_line.svg'),'-dsvg',print_res) 
    print(fig,fullfile(output_result_path,'/p_hd_x1_line.png'),'-dpng',print_res)
    savefig(fig,fullfile(output_result_path,'/p_hd_x1_line.fig'))
end
clear fig; clear ax;
clear x2_0; clear x2_0_i;
% Surface:
fig = figure('Units','centimeters','Position',[12.1 01 figuresize]);
surf(h.x1*1e6,h.x2*1e6,sol.p_hd'*1e-6)
shading interp
material dull
colormap(colmap_p)
camlight
ax = gca;
ax.FontSize = sizeofticksfonts;
xlim(x1_range)
xticks(x1_ticks)
ylim(x2_range)
yticks(x2_ticks)
zlim(p_range)
zticks(p_ticks)
xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
zlabel('$p_{hd} \mathrm{[MPa]}$','fontsize',sizeoffonts)
if flag_save_plots 
    print(fig,fullfile(output_result_path,'/p_hd.svg'),'-dsvg',print_res) 
    print(fig,fullfile(output_result_path,'/p_hd.png'),'-dpng',print_res)
    savefig(fig,fullfile(output_result_path,'/p_hd.fig'))
end
clear fig; clear ax;

% Fluid shear stress in positive x1-direction acting from the fluid upon the upper surface 
fig = figure('Units','centimeters','Position',[18.1 01 figuresize]);
surf(h.x1*1e6,h.x2*1e6,sol.tau_hd_up'*1e-6)
shading interp
material dull
colormap(colmap_p)
camlight
ax = gca;
ax.FontSize = sizeofticksfonts;
xlim(x1_range)
xticks(x1_ticks)
ylim(x2_range)
yticks(x2_ticks)
zlim(flip(-tau_range))
zticks(flip(-tau_ticks))
xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
zlabel('$\tau_{hd,up} \mathrm{[MPa]}$','fontsize',sizeoffonts)
if flag_save_plots 
    print(fig,fullfile(output_result_path,'/tau_hd_up.svg'),'-dsvg',print_res) 
    print(fig,fullfile(output_result_path,'/tau_hd_up.png'),'-dpng',print_res)
    savefig(fig,fullfile(output_result_path,'/tau_hd_up.fig'))
end
clear fig; clear ax;

% Fluid shear stress in positive x1-direction acting from the fluid upon the lower surface 
fig = figure('Units','centimeters','Position',[18.1 01 figuresize]);
surf(h.x1*1e6,h.x2*1e6,sol.tau_hd_low'*1e-6)
shading interp
material dull
colormap(colmap_p)
camlight
ax = gca;
ax.FontSize = sizeofticksfonts;
xlim(x1_range)
xticks(x1_ticks)
ylim(x2_range)
yticks(x2_ticks)
zlim(tau_range)
zticks(tau_ticks)
xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
zlabel('$\tau_{hd,low} \mathrm{[MPa]}$','fontsize',sizeoffonts)
if flag_save_plots 
    print(fig,fullfile(output_result_path,'/tau_hd_low.svg'),'-dsvg',print_res) 
    print(fig,fullfile(output_result_path,'/tau_hd_low.png'),'-dpng',print_res)
    savefig(fig,fullfile(output_result_path,'/tau_hd_low.fig'))
end
clear fig; clear ax;

% Cavity fraction:
fig = figure('Units','centimeters','Position',[18.1 01 figuresize]);
surf(h.x1*1e6,h.x2*1e6,sol.t')
shading interp
material dull
colormap(colmap_p)
camlight
ax = gca;
ax.FontSize = sizeofticksfonts;
xlim(x1_range)
xticks(x1_ticks)
ylim(x2_range)
yticks(x2_ticks)
zlim(t_range)
zticks(t_ticks)
xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
zlabel('$\theta \mathrm{[-]}$','fontsize',sizeoffonts)
if flag_save_plots 
    print(fig,fullfile(output_result_path,'/t.svg'),'-dsvg',print_res) 
    print(fig,fullfile(output_result_path,'/t.png'),'-dpng',print_res)
    savefig(fig,fullfile(output_result_path,'/t.fig'))
end
clear fig; clear ax;

% Gap height:
% Line along x1:
[x2_0,x2_0_i] = min(abs(h.x2));
fig = figure('Units','centimeters','Position',[00 11 figuresize]);
plot(h.x1*1e6,h.h_ma_ad(:,x2_0_i)*1e9,'LineWidth',widthlines,'Color',KIT_colorlist{1})
grid on
ax = gca;
ax.FontSize = sizeofticksfonts;
xlim(x1_range)
xticks(x1_ticks)
ylim(h_range)
yticks(h_ticks)
xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
ylabel('$h_{ma,ad} \mathrm{[nm]}$','fontsize',sizeoffonts);
title("$x_2=$ " + x2_0*1e6 + "$\mathrm{\mu m}$",'fontsize',sizeoffonts)
if flag_save_plots 
    print(fig,fullfile(output_result_path,'/Gap_height_ma_adj_x1_line.svg'),'-dsvg',print_res) 
    print(fig,fullfile(output_result_path,'/Gap_height_ma_adj_x1_line.png'),'-dpng',print_res)
    savefig(fig,fullfile(output_result_path,'/Gap_height_ma_adj_x1_line.fig'))
end
clear fig; clear ax;
clear x2_0; clear x2_0_i;
% Line along x2:
[x1_0,x1_0_i] = min(abs(h.x1));
fig = figure('Units','centimeters','Position',[00 11 figuresize]);
plot(h.x2*1e6,h.h_ma_ad(x1_0_i,:)*1e9,'LineWidth',widthlines,'Color',KIT_colorlist{1})
grid on
ax = gca;
ax.FontSize = sizeofticksfonts;
xlim(x2_range)
xticks(x2_ticks)
ylim(h_range)
yticks(h_ticks)
xlabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
ylabel('$h_{ma,ad} \mathrm{[nm]}$','fontsize',sizeoffonts);
title("$x_1=$ " + x1_0*1e6 + "$\mathrm{\mu m}$",'fontsize',sizeoffonts)
if flag_save_plots 
    print(fig,fullfile(output_result_path,'/Gap_height_ma_adj_x2_line.svg'),'-dsvg',print_res) 
    print(fig,fullfile(output_result_path,'/Gap_height_ma_adj_x2_line.png'),'-dpng',print_res)
    savefig(fig,fullfile(output_result_path,'/Gap_height_ma_adj_x2_line.fig'))
end
clear fig; clear ax;
clear x1_0; clear x1_0_i;
% Contour:
[x1_matr, x2_matr] = ndgrid(h.x1,h.x2);
fig = figure('Units','centimeters','Position',[00 11 figuresize]);
[~,l] = contourf(x1_matr'*1e6,x2_matr'*1e6,h.h_ma_ad'*1e9,h_cont_ticks);
l.LineWidth = contour_line_width;
axis equal 
ax = gca;
ax.FontSize = sizeofticksfonts;
c = colorbar('TickLabelInterpreter','latex');
c.Ticks = h_ticks;
c.FontSize = sizeofticksfonts;
c.Label.Interpreter = 'latex';
c.Label.String = '$h_{ma,ad} \mathrm{[nm]}$';
c.Label.FontSize = sizeoffonts;
colormap(flipud(colmap_h))
xlim(x1_range)
xticks(x1_ticks)
ylim(x2_range)
yticks(x2_ticks)
caxis(h_cont_range)
xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
if flag_save_plots 
    print(fig,fullfile(output_result_path,'/Gap_height_ma_adj_cont.svg'),'-dsvg',print_res) 
    print(fig,fullfile(output_result_path,'/Gap_height_ma_adj_cont.png'),'-dpng',print_res)
    savefig(fig,fullfile(output_result_path,'/Gap_height_ma_adj_cont.fig'))
end
clear fig; clear x1_matr; clear x2_matr; clear c; clear ax; clear l;
% Surface:
fig = figure('Units','centimeters','Position',[00 11 figuresize]);
surf(h.x1*1e6,h.x2*1e6,h.h_ma_ad'*1e9)
fig.CurrentAxes.ZDir = 'Reverse';
shading interp
material dull
colormap(flipud(colmap_h))
camlight
ax = gca;
ax.FontSize = sizeofticksfonts;
xlim(x1_range)
xticks(x1_ticks)
ylim(x2_range)
yticks(x2_ticks)
zlim(h_range)
caxis(h_range)
xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
zlabel('$h_{ma,ad} \mathrm{[nm]}$','fontsize',sizeoffonts);
if flag_save_plots 
    print(fig,fullfile(output_result_path,'/Gap_height_ma_adj.svg'),'-dsvg',print_res) 
    print(fig,fullfile(output_result_path,'/Gap_height_ma_adj.png'),'-dpng',print_res)
    savefig(fig,fullfile(output_result_path,'/Gap_height_ma_adj.fig'))
end
clear fig; clear ax;

% Total pressure:
fig = figure('Units','centimeters','Position',[00 01 figuresize]);
surf(h.x1*1e6,h.x2*1e6,sol.p_tot'*1e-6)
shading interp
material dull
colormap(colmap_p)
camlight
ax = gca;
ax.FontSize = sizeofticksfonts;
xlim(x1_range)
xticks(x1_ticks)
ylim(x2_range)
yticks(x2_ticks)
zlim(p_range)
zticks(p_ticks)
xlabel('$x_1 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
ylabel('$x_2 \mathrm{[\mu m]}$','fontsize',sizeoffonts);
zlabel('$p_{tot} \mathrm{[MPa]}$','fontsize',sizeoffonts)
if flag_save_plots 
    print(fig,fullfile(output_result_path,'/p_tot.svg'),'-dsvg',print_res) 
    print(fig,fullfile(output_result_path,'/p_tot.png'),'-dpng',print_res)
    savefig(fig,fullfile(output_result_path,'/p_tot.fig'))
end
clear fig; clear ax;

set(groot,'defaulttextinterpreter','remove');  
set(groot,'defaultAxesTickLabelInterpreter','remove');  
set(groot,'defaultLegendInterpreter','remove');

%% Write table:
Input_Property = {'p_amb'; 'p_cav'; 'rho_0'; 'mu_0'; 'alpha_mu'; ...
    'nu_up'; 'nu_low'; 'E_up'; 'E_low'; 'E_dash'; ...
    'W'; 'u_up'; 'u_low'};
Unit = {'[Pa]'; '[Pa]'; '[kg/m^3]'; '[Pas]'; '[-]'; ...
    '[-]'; '[-]'; '[-]'; '[Pa]'; '[Pa]'; ...
    '[N]'; '[m/s]'; '[m/s]'};
Value    = [fld.p_amb; fld.p_cav; fld.rho_0; fld.mu_0; fld.alpha; ...
    sld.nu_up; sld.nu_low; sld.E_up; sld.E_low; sld.E_dash; ...
    opc.W; opc.u_up(i_OC); opc.u_low];
Table_Input_Properties = table(Input_Property, Unit, Value);
clear Input_Property; clear Unit; clear Value;
if flag_save_plots 
    writetable(Table_Input_Properties,fullfile(output_result_table_path,'/Input_properties.csv'),'Delimiter',',','QuoteStrings',true)
end