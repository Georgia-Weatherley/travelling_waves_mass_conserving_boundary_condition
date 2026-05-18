
clear; clc; close all;
% This script reproduces Figure 7 from the manuscript. It plots an invading
% travelling wave and a receding travelling wave 

% It calls the following functions:  
% solve_phase_plane.m : solves the ODE system 
% Solve_PDE_func.m : identifies the front density u*

% It also loads the .mat file 
% wave_speed_curve.mat : data for the wave speed curve 

%% Plot styles 
figure
tiledlayout(2,2);

Traj_LW   = 5;
Spiral_LW = 3;
M_size    = 12;
Arrow_LW  = 1.5;
Axis_FS   = 30;
Axis_LW   = 1.5;
BC_LW     = 3;

BC_col   = 1/255*[39,73,84];
traj_col = 1/255*[136 65 157];

subtle_red  = 1/255*[254,196,79];
yellow_col  = 1/255*[254,196,79];
Light_purple = 1/255*[158,188,218];
Light_green  = 1/255*[103,169,207];
Grey_col     = 1/255*[115,115,115];

% Load wave speed curve
load("wave_speed_curve.mat")

%% Invading case
a = 0.5;
b = 1;

BC_func = @(u) a*u + b;

BC_vals   = BC_func(u_star_PP);
difference = BC_vals - c_array;

idx = find(difference(1:end-1).*difference(2:end) < 0);
idx = idx(1);

x_intersect = interp1(difference(idx:idx+1),u_star_PP(idx:idx+1),0);

y_intersect = BC_func(x_intersect);

c = y_intersect;

sol_ode = solve_phase_plane(c);

t_ode = sol_ode.x;
U_ode = sol_ode.y(1,:);

u_shift = interp1(U_ode, t_ode, x_intersect);

indx = find(t_ode > u_shift);
indx = indx(1);

x_val_ODE = t_ode(indx-1);

[t,sol,rho] = Solve_PDE_func( ...
    a,b,25,10000,0.5,0,200,401,[],0);

% Plot PDE sol
nexttile(1)
hold on

t_vals = [2 10 20 30 40 50 60];

time_index = zeros(size(t_vals));

for kk = 1:length(t_vals)
    time_index(kk) = find(t == t_vals(kk));
end

LL = sol(:,end);

plot([0 200],[x_intersect x_intersect],':','LineWidth',BC_LW,'Color','k')

plot([rho'*LL(1) LL(1)],[sol(1,1:end-1) 0],'-','LineWidth',Traj_LW,'Color',traj_col)

for kk = 1:length(t_vals)
    plot(rho.*LL(time_index(kk)),sol(time_index(kk),1:end-1),'-','LineWidth',Traj_LW,'Color',traj_col)

end

% Overlay travelling wave from phase plane 

x_val_last_PDE = rho(end).*LL(time_index(end));

spatial_shift = x_val_last_PDE - x_val_ODE;

t_fill = linspace(0,t_ode(1)+spatial_shift,100);

plot([t_fill t_ode(1:indx-1)+spatial_shift],[ones(length(t_fill),1); U_ode(1:indx-1)'], ...
    ':','Color',yellow_col,'LineWidth',Traj_LW)

xlabel('$x$','Interpreter','latex')
ylabel('$u(x,t)$','Interpreter','latex')

set(gca,'FontSize',Axis_FS, ...
    'TickLabelInterpreter','latex')

xlim([0 100])

grid off, axis on, box on

% Plot the boundary speed 

nexttile(2)
hold on

dL_dt = a*sol(:,end-1) + b;

plot(t(1:end),dL_dt,'-','LineWidth',Traj_LW,'Color',Grey_col)

plot([0 200],[y_intersect y_intersect],':','LineWidth',BC_LW,'Color','r')

for kk = 1:length(t_vals)

    plot(t(time_index(kk)),dL_dt(time_index(kk)),'diamond','MarkerSize',M_size, ...
         'MarkerEdgeColor',traj_col,'MarkerFaceColor',traj_col)

end

xlim([0 65])

xlabel('$t$','Interpreter','latex')
ylabel('$dL/dt$','Interpreter','latex')

set(gca,'FontSize',Axis_FS,'TickLabelInterpreter','latex')

grid off, axis on, box on

%% Receding case 
a = -0.3;
b = -0.2;

BC_func = @(u) a*u + b;

BC_vals = BC_func(u_star_PP);

difference = BC_vals - c_array;

idx = find(difference(1:end-1).*difference(2:end) < 0);

x_intersect = interp1(difference(idx:idx+1),u_star_PP(idx:idx+1),0);

y_intersect = BC_func(x_intersect);

c = y_intersect;

sol_ode = solve_phase_plane(c);

t_ode = sol_ode.x;
U_ode = sol_ode.y(1,:);

u_shift = interp1(U_ode, t_ode, x_intersect);

indx = find(t_ode > u_shift);
indx = indx(1);

x_val_ODE = t_ode(indx-1);

traj_col = 1/255*[1 108 89];

% Solve PDE
[t,sol,rho] = Solve_PDE_func(a,b,100,10000,0.5,0,200,401,[],0);

% Plot PDE sol
nexttile(3)
hold on

t_vals = [2 10 20 30 40 50 60];

time_index = zeros(size(t_vals));

for kk = 1:length(t_vals)
    time_index(kk) = find(t == t_vals(kk));
end

LL = sol(:,end);

plot([0 1000],[x_intersect x_intersect],':','LineWidth',BC_LW,'Color','k')

plot([rho'*LL(1) LL(1)],[sol(1,1:end-1) 0],'-','LineWidth',Traj_LW,'Color',traj_col)

for kk = 1:length(t_vals)

    plot(rho.*LL(time_index(kk)),sol(time_index(kk),1:end-1),'-','LineWidth',Traj_LW,'Color',traj_col)

end

% Overlay travelling wave from phase plane 

x_val_last_PDE = rho(end).*LL(time_index(end));

spatial_shift = x_val_last_PDE - x_val_ODE;

t_fill = linspace(0,t_ode(1)+spatial_shift,100);

plot([t_fill, t_ode(1:indx)+spatial_shift],[ones(length(t_fill),1); U_ode(1:indx)'],':', ...
    'Color',yellow_col, 'LineWidth',Traj_LW)

xlabel('$x$','Interpreter','latex')
ylabel('$u(x,t)$','Interpreter','latex')

set(gca,'FontSize',Axis_FS, ...
    'TickLabelInterpreter','latex')

xlim([0 100])

grid off, axis on, box on

% Plot the boundary speed 

nexttile(4)
hold on

dL_dt = a*sol(:,end-1) + b;

plot(t(1:end),dL_dt,'-','LineWidth',Traj_LW,'Color',Grey_col)

plot([0 200],[y_intersect y_intersect],':','LineWidth',BC_LW,'Color','r')

for kk = 1:length(t_vals)

    plot(t(time_index(kk)),dL_dt(time_index(kk)),'diamond','MarkerSize',M_size, ...
        'MarkerEdgeColor',traj_col,'MarkerFaceColor',traj_col)

end

xlim([0 65])

xlabel('$t$','Interpreter','latex')
ylabel('$dL/dt$','Interpreter','latex')

set(gca,'FontSize',Axis_FS,'TickLabelInterpreter','latex')

grid off, axis on, box on