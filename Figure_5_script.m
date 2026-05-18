clear; clc; close all;
% This script reproduces Figure 5 from the manuscript. It shows
% intersections of a linear boundary velocity function with the wave speed
% curve. 

% It loads the .mat file 
% wave_speed_curve.mat : data for the wave speed curve 

load("wave_speed_curve.mat")

%% Plot styles
BC_LW = 4; 
Traj_LW = 5;
Axis_LW   = 1.5;
Axis_FS   = 30;
traj_col = [0.2431    0.2980    0.3686];

map = 1/255*[77,0,75; 129,15,124; 136,65,157; 140,107,177; 140,150,198; 158,188,218; 191,211,230; 166,189,219; 103,169,207; 54,144,192; 2,129,138; 1,108,89; 1,70,54];
grey_col = 1/255*[220 220 220]; 
BC_col = 1/255*[223,101,176];
BC_col = 1/255*[201,148,199];
BC_col = 1/255*[247,104,161];
BC_col = 1/255*[190 190 190]; 
U_array = linspace(-1.7, 8 ,200);
BC_func = @(u,a,b) a*u + b;

figure, tiledlayout(1,3)

%% Subplot 1: one invading 
nexttile(1), hold on 
alpha = 0.5;
beta = 1;

purple = 1/255*[216,199,224];
purple = 1/255*[136 65 157];
green = 1/255*[185,210,205];
green = 1/255*[1 108 89];
M_size = 20;

S = 0.5; L = 1.5;
run  = L / sqrt(1 + S^2); rise = S * run;
theta1 = atan(rise/run);
x1 = 2; y1 = 1.4;

x = [x1 (x1-run) x1 x1]; 
y = [y1 (y1-rise) (y1-rise) y1];

fill(x,y,grey_col,EdgeColor="none")
plot([-5 50],[0 0],'Linewidth',Axis_LW,'Color','k','HandleVisibility','off'); plot([-5 50],[0 0],'Linewidth',Axis_LW,'Color','k','HandleVisibility','off');
plot([0 0],[-50 50],'Linewidth',Axis_LW,'Color','k','HandleVisibility','off'); plot([-5 50],[0 0],'Linewidth',Axis_LW,'Color','k','HandleVisibility','off');
plot(u_star_PP,c_array,'-','Color',traj_col,'LineWidth',Traj_LW,'HandleVisibility','off');
plot(U_array,BC_func(U_array,alpha,beta),'LineWidth',BC_LW,'Color',BC_col)
plot(0.1698,1.0849,'o','MarkerFaceColor',purple,'MarkerEdgeColor',purple,'MarkerSize',M_size)
h = text(-0.485, 1.1, '$\beta = 1$','Interpreter','latex','FontSize',24);
h = text(((min(x) + max(x))/2) -.25, ((min(y) + max(y))/2), '$\alpha = 0.5$','Interpreter','latex','FontSize',24);
set(h,'Rotation',rad2deg(theta1));
xlim([-0.5 2.5]), ylim([-1 2.2])
axis square 
xlabel('$U^{*}$',Interpreter='latex'), ylabel('c',Interpreter='latex')
set(gca,'FontSize',Axis_FS, 'TickLabelInterpreter', 'latex')
box on 

%% Subplot 2: one receding 
nexttile(2), hold on  
alpha = -0.3;
beta = -0.2;

S = -0.3;
L = 1.5;
run  = L / sqrt(1 + S^2); rise = S * run;
theta2 = atan(rise/run);
x1 = 0.8; y1 = 1.4;

x2 = [x1 x1 (x1+run) x1];
y2 = [y1 (y1+rise) (y1+rise) y1];
fill(x2,y2,grey_col,EdgeColor="none")

plot([-5 50],[0 0],'Linewidth',Axis_LW,'Color','k','HandleVisibility','off'); plot([-5 50],[0 0],'Linewidth',Axis_LW,'Color','k','HandleVisibility','off');
plot([0 0],[-50 50],'Linewidth',Axis_LW,'Color','k','HandleVisibility','off'); plot([-5 50],[0 0],'Linewidth',Axis_LW,'Color','k','HandleVisibility','off');
plot(u_star_PP,c_array,'-','Color',traj_col,'LineWidth',Traj_LW,'HandleVisibility','off');
plot(U_array,BC_func(U_array,alpha,beta),'LineWidth',BC_LW,'Color',BC_col)
plot(1.8246,-0.7474,'o','MarkerFaceColor',green,'MarkerEdgeColor',green,'MarkerSize',M_size)
h = text(0.16, -0.15, '$\beta = -0.2$','Interpreter','latex','FontSize',24);
h = text(((min(x) + max(x))/2) -.25, ((min(y) + max(y))/2)+0.35, '$\alpha = -0.3$','Interpreter','latex','FontSize',24);
set(h,'Rotation',rad2deg(theta2));
xlim([-0.5 2.5]), ylim([-1 2.2])
axis square
xlabel('$U^{*}$',Interpreter='latex'), ylabel('c',Interpreter='latex')
set(gca,'FontSize',Axis_FS, 'TickLabelInterpreter', 'latex')
box on 

%% Subplot 3: one invading, one receding 
nexttile(3), hold on 
alpha = -1.4;
beta = 1.5;

x2 = 2; y2 = 1.4;
S = -1.4; L = 1.5;
run  = L / sqrt(1 + S^2); rise = S * run;
theta3 = atan(rise/run);

x = [x2 (x2-run) x2 x2]; 
y = [y2 y2 (y2+rise) y2];

fill(x,y,grey_col,EdgeColor="none")

plot([-5 50],[0 0],'Linewidth',Axis_LW,'Color','k','HandleVisibility','off'); plot([-5 50],[0 0],'Linewidth',Axis_LW,'Color','k','HandleVisibility','off');
plot([0 0],[-50 50],'Linewidth',Axis_LW,'Color','k','HandleVisibility','off'); plot([-5 50],[0 0],'Linewidth',Axis_LW,'Color','k','HandleVisibility','off');
plot(u_star_PP,c_array,'-','Color',traj_col,'LineWidth',Traj_LW,'HandleVisibility','off');
plot(U_array,BC_func(U_array,alpha,beta),'LineWidth',BC_LW,'Color',BC_col)
plot(1.2305, -0.2227,'o','MarkerFaceColor',green,'MarkerEdgeColor',green,'MarkerSize',M_size) 
plot(0.0446, 1.4376,'o','MarkerFaceColor',purple,'MarkerEdgeColor',purple,'MarkerSize',M_size)
h = text(0.2, 1.5, '$\beta = 1.5$','Interpreter','latex','FontSize',24);
h = text(((min(x) + max(x))/2)-0.3, ((min(y) + max(y))/2)+0.2, '$\alpha = -1.4$','Interpreter','latex','FontSize',24);
set(h,'Rotation',rad2deg(theta3));
xlim([-0.5 2.5]), ylim([-1 2.2])
axis square
xlabel('$U^{*}$',Interpreter='latex'), ylabel('c',Interpreter='latex')
set(gca,'FontSize',Axis_FS, 'TickLabelInterpreter', 'latex')
box on 
