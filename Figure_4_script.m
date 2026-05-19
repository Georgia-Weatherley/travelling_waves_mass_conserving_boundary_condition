clear; clc; close all;
% This script reproduces Figure 4 from the manuscript. 
% It calls the functions: 
% solve_phase_plane.m 
% Find_ustar_func.m
% Load_pert_sols.m

% It also loads the wave speed curve (u* vs c) from:
% "wave_speed_curve.mat"

clear; close all; clc;

load("wave_speed_curve.mat")

pos_traj_col = 1/255*[136 65 157];

pos_pert_col  = 1/255*[238 109 113];
pos_pert_col2 = 1/255*[203 24 29];

BC_col      = 1/255*[39 73 84];
green_neg4  = 1/255*[1 108 89];

traj_col = [0.2431 0.2980 0.3686];

peach      = 1/255*[247 197 159];
dark_peach = 1/255*[208 99 17];

Traj_LW = 4;
BC_LW   = 3;
Axis_LW = 1.5;
Axis_FS = 30;
M_size = 20;

%% Plot styles 
all_cols = [0.5059    0.0588    0.4863;
    0.5333    0.2549    0.6157;
    0.5490    0.4196    0.6941;
    0.5490    0.5882    0.7765;
    0.6510    0.7412    0.8588;
    0.4039    0.6627    0.8118;
    0.2118    0.5647    0.7529;
    0.0078    0.5059    0.5412;
    0.0039    0.4235    0.3490;
    0.0039    0.2745    0.2118];


%% Asymptotic approximations
[W_pert_neg_c, W_pert_small_c,c_equ_neg_c, c_equ_small_c] = load_pert_sols;

zmax = 40;
npts = 4000;

%%  FIGURE 1 : Negative wave speeds
figure
tiledlayout(1,3)

%%  c = -3

nexttile(1)
hold on

c = -3;
w_func = @(u) -c*u;

sol = solve_phase_plane(c);
z = linspace(0,zmax,npts);
Y = interp1(sol.x.',sol.y.',z,'spline').';

U = Y(1,:)';
V = Y(2,:)';

[u_intercept,~] = Find_ustar_func(c,U,V,w_func,[]);

spec_point_neg1 = [c u_intercept];

U_array = 1:0.01:8;

% Axes
plot([0 0],[-1 120],'LineWidth',Axis_LW,'Color','k')
plot([-0.5 20],[0 0],'LineWidth',Axis_LW,'Color','k')

% Numerical trajectory
plot(U,V,'-','LineWidth',Traj_LW,'Color',all_cols(8,:))

% Asymptotic approximation
plot(U_array,W_pert_neg_c(U_array,c),':','LineWidth',Traj_LW,'Color',peach)

% Boundary condition line
BC_array = linspace(0,u_intercept,200);

plot(BC_array,w_func(BC_array),'-','LineWidth',BC_LW,'Color',BC_col)

% Intersection point
plot(u_intercept,w_func(u_intercept),'o','MarkerSize',M_size,'MarkerFaceColor',all_cols(8,:), ...
    'MarkerEdgeColor',all_cols(8,:))

% Equil
plot(1,0,'ko','MarkerSize',M_size,'MarkerFaceColor','k')

plot(0,0,'ko','MarkerSize',M_size,'MarkerFaceColor','k')

xlabel('$U(z)$','Interpreter','latex')
ylabel('$W(z)$','Interpreter','latex')

set(gca,'FontSize',Axis_FS,'TickLabelInterpreter','latex')

xlim([-0.2 6])
ylim([-0.5 20])

axis square, box on, grid off

%%  c = -5
nexttile(2)
hold on

c = -5;
w_func = @(u) -c*u;

sol = solve_phase_plane(c);

z = linspace(0,zmax,npts);

Y = interp1(sol.x.',sol.y.',z,'spline').';

U = Y(1,:)';
V = Y(2,:)';

[u_intercept,~] = Find_ustar_func(c,U,V,w_func,[]);

spec_point_neg2 = [c u_intercept];

% Asymptotic curve
U_array = 1:0.01:16;

% Axes
plot([0 0],[-1 120],'LineWidth',Axis_LW,'Color','k')

plot([-0.5 20],[0 0],'LineWidth',Axis_LW,'Color','k')

% Numerical trajectory
plot(U,V,'-','LineWidth',Traj_LW,'Color',all_cols(10,:))

% Asymptotic approximation
plot(U_array,W_pert_neg_c(U_array,c),':','LineWidth',Traj_LW,'Color',peach)

% Boundary condition line
BC_array = linspace(0,u_intercept,200);

plot(BC_array,w_func(BC_array),'-','LineWidth',BC_LW,'Color',BC_col)

% Intersection point
plot(u_intercept,w_func(u_intercept),'o','MarkerSize',M_size, ...
    'MarkerFaceColor',green_neg4,'MarkerEdgeColor',all_cols(10,:))

% Equilibria
plot(1,0,'ko','MarkerSize',M_size,'MarkerFaceColor','k')

plot(0,0,'ko','MarkerSize',M_size, ...
    'MarkerFaceColor','k')

xlabel('$U(z)$','Interpreter','latex')
ylabel('$W(z)$','Interpreter','latex')

set(gca,...
    'FontSize',Axis_FS,...
    'TickLabelInterpreter','latex')

xlim([-0.3 10])
ylim([-1.2 60])

axis square
box on
grid off

%% Negative c asymptotic comparison

nexttile(3)
hold on

plot([-5 50],[0 0],'LineWidth',Axis_LW,'Color','k')

plot([0 0],[-50 50],'LineWidth',Axis_LW,'Color','k')

% Numerical branch
plot(u_star_PP,c_array,'-','LineWidth',Traj_LW,'Color',traj_col)

% Asymptotic branch
plot(u_star_PP,c_equ_neg_c(u_star_PP),':','LineWidth',Traj_LW,'Color',dark_peach)

% Highlight selected points
plot(spec_point_neg1(2),spec_point_neg1(1),'o','MarkerSize',M_size, ...
    'MarkerFaceColor',all_cols(8,:),'MarkerEdgeColor',all_cols(8,:))

plot(spec_point_neg2(2),spec_point_neg2(1),'o','MarkerSize',M_size, ...
    'MarkerFaceColor',all_cols(10,:),'MarkerEdgeColor',all_cols(10,:))

xlabel('$U^{*}$','Interpreter','latex')
ylabel('$c$','Interpreter','latex')

set(gca,...
    'FontSize',Axis_FS,...
    'TickLabelInterpreter','latex')

xlim([-0.5 14.5])
ylim([-10 2])

axis square, box on, grid off

%% FIGURE 2 : Small positive wave speeds

figure
tiledlayout(1,3)

%%  c = 0.5

nexttile(1)
hold on

c = 0.5;
w_func = @(u) -c*u;

sol = solve_phase_plane(c);

z = linspace(0,zmax,npts);

Y = interp1(sol.x.',sol.y.',z,'spline').';

U = Y(1,:)';
V = Y(2,:)';

[u_intercept,~] = Find_ustar_func(c,U,V,w_func,[]);
[~,u_below2]    = Find_ustar_func(c,U,V,[],0);

spec_point_pos1 = [c u_intercept];

% Asymptotic curve
U_array = 0:0.01:1;

% Axes
plot([0 0],[-1 25],'LineWidth',Axis_LW,'Color','k')

plot([-0.5 20],[0 0],'LineWidth',Axis_LW,'Color','k')

% Numerical trajectory
plot(U(1:u_below2),V(1:u_below2),'-','LineWidth',Traj_LW,'Color',all_cols(4,:))

% Asymptotic approximation
plot(U_array,W_pert_small_c(U_array,c),':','LineWidth',Traj_LW,'Color',pos_pert_col2)

% Boundary condition
BC_array = linspace(0,u_intercept,200);

plot(BC_array,w_func(BC_array),'-','LineWidth',BC_LW,'Color',BC_col)

% Intersection point
plot(u_intercept,w_func(u_intercept),'o','MarkerSize',M_size, ...
    'MarkerFaceColor',all_cols(4,:),'MarkerEdgeColor',all_cols(4,:))

% Equilibria
plot(1,0,'ko','MarkerSize',M_size,'MarkerFaceColor','k')

plot(0,0,'ko','MarkerSize',M_size,'MarkerFaceColor','k')

xlabel('$U(z)$','Interpreter','latex')
ylabel('$W(z)$','Interpreter','latex')

set(gca,'FontSize',Axis_FS,'TickLabelInterpreter','latex')

xlim([-0.05 1.1])
ylim([-0.5 0.05])

axis square, box on, grid off

%%  c = 0.25

nexttile(2)
hold on

c = 0.25;
w_func = @(u) -c*u;

sol = solve_phase_plane(c);

z = linspace(0,zmax,npts);

Y = interp1(sol.x.',sol.y.',z,'spline').';

U = Y(1,:)';
V = Y(2,:)';

[u_intercept,~] = Find_ustar_func(c,U,V,w_func,[]);
[~,u_below2]    = Find_ustar_func(c,U,V,[],0);

spec_point_pos2 = [c u_intercept];

% Asymptotic curve
U_array = 0:0.01:1;

% Axes
plot([0 0],[-1 25],'LineWidth',Axis_LW,'Color','k')

plot([-0.5 20],[0 0],'LineWidth',Axis_LW,'Color','k')

% Numerical trajectory
plot(U(1:u_below2),V(1:u_below2),'-','LineWidth',Traj_LW,'Color',all_cols(3,:))

% Asymptotic approximation
plot(U_array,W_pert_small_c(U_array,c),':','LineWidth',Traj_LW,'Color',pos_pert_col2)

% Boundary condition
BC_array = linspace(0,u_intercept,200);

plot(BC_array,w_func(BC_array),'-','LineWidth',BC_LW,'Color',BC_col)

% Intersection point
plot(u_intercept,w_func(u_intercept),'o','MarkerSize',M_size,'MarkerFaceColor',all_cols(3,:), ...
    'MarkerEdgeColor',all_cols(3,:))

% Equilibria
plot(1,0,'ko','MarkerSize',M_size,'MarkerFaceColor','k')

plot(0,0,'ko','MarkerSize',M_size,'MarkerFaceColor','k')

xlabel('$U(z)$','Interpreter','latex')
ylabel('$W(z)$','Interpreter','latex')

set(gca,'FontSize',Axis_FS,'TickLabelInterpreter','latex')

xlim([-0.05 1.1])
ylim([-0.5 0.05])

axis square, box on, grid off

%%  Small positive c asymptotic comparison

nexttile(3)
hold on

plot([0 0],[-50 50],'LineWidth',Axis_LW,'Color','k')

plot([-5 50],[0 0],'LineWidth',Axis_LW,'Color','k')

% Numerical branch
plot(u_star_PP,c_array,'-','LineWidth',Traj_LW,'Color',traj_col)

% Asymptotic branch
plot(u_star_PP,c_equ_small_c(u_star_PP),':','LineWidth',Traj_LW,'Color',pos_pert_col)

% Highlight selected points
plot(spec_point_pos1(2),spec_point_pos1(1),'o','MarkerSize',M_size, ...
    'MarkerFaceColor',all_cols(4,:),'MarkerEdgeColor',all_cols(4,:))

plot(spec_point_pos2(2),spec_point_pos2(1),'o','MarkerSize',M_size, ...
    'MarkerFaceColor',all_cols(3,:),'MarkerEdgeColor',all_cols(3,:))

xlabel('$U^{*}$','Interpreter','latex')
ylabel('$c$','Interpreter','latex')

set(gca,'FontSize',Axis_FS,'TickLabelInterpreter','latex')

xlim([-0.1 2])
ylim([-1 2])

axis square, box on, grid off