clear; clc; close all;
% This script reproduces Figure 2 from the manuscript. 

% It calls the following functions:  
% solve_phase_plane.m : solves the ODE system 
% Find_ustar_func.m : identifies the front density u*

% It contains the following functions:  
% plot_profile.m : plots the density profile of a travelling wave solution 
% plot_case.m : plots the phase portrait of a travelling wave solution 

%%  Plot style

Traj_LW   = 3;
Spiral_LW = 2;
M_size    = 18;
Arrow_LW  = 1.5;
Axis_FS   = 30;
Axis_LW   = 1.5;
BC_LW     = 2;

BC_col = 1/255*[39,73,84];
quiv_col = 1/255*[149 161 178];

%%  Define and plot positive wave speeds 

% Cols for positive values of c 
traj_col = 1/255*[136 65 157];
Our_BC_col  = 1/255*[159 82 183];
Maud_BC_col = 1/255*[253 141 60];
FS_BC_col   = 1/255*[217 72 1];

figure
tl = tiledlayout(2,2);

% Plot phase portrait 
plot_case(0.25,1)
plot_case(1,3)

% Plot density profile 
plot_profile(0.25,2)
plot_profile(1,4)

%%  Define and plot negative wave speeds 

% Cols for negative values of  c
traj_col = 1/255*[1 108 89];
Our_BC_col  = 1/255*[1 141 116];

figure
tl = tiledlayout(2,2);

% Plot phase portrait 
plot_case(-1.5,1)
plot_case(-3,3)

% Plot density profile 
plot_profile(-1.5,2)
plot_profile(-3,4)

%%  Phase plane plotting function
function plot_case(c,tile_num)

    % Plot Styles
    Traj_LW  = evalin('base','Traj_LW');
    Arrow_LW = evalin('base','Arrow_LW');
    Axis_FS  = evalin('base','Axis_FS');
    Axis_LW  = evalin('base','Axis_LW');
    BC_LW    = evalin('base','BC_LW');
    M_size   = evalin('base','M_size');

    BC_col       = evalin('base','BC_col');
    quiv_col     = evalin('base','quiv_col');
    traj_col     = evalin('base','traj_col');
    Our_BC_col   = evalin('base','Our_BC_col');
    FS_BC_col    = evalin('base','FS_BC_col');
    Maud_BC_col  = evalin('base','Maud_BC_col');

    % Solve ODE system
    sol = solve_phase_plane(c);

    U = sol.y(1,:)';
    V = sol.y(2,:)';

    % Boundary function
    w_func = @(u) -c*u;

    % Intersections
    [u_intercept,u_below] = Find_ustar_func(c,U,V,w_func,[]);

    if c > 0
        [v_val_FS,~]   = Find_ustar_func(c,U,V,[],0);
        [v_val_Maud,~] = Find_ustar_func(c,U,V,[],0.5);
    end

    nexttile(tile_num)
    hold on

    % Vector field
    if c > 0

        [u,v] = meshgrid(linspace(-0.5,1.5,10),linspace(-0.5,0.5,10));

    else

        [u,v] = meshgrid(linspace(-1,8,10),linspace(-1,25,10));

    end

    du = v;
    dv = -c.*v - u.*(1-u);

    quiver(u,v,du,dv,1,'LineWidth',Arrow_LW,'Color',[quiv_col 0.2])

    % Axes
    if c > 0

        plot([0 0],[-0.5 2],'k','LineWidth',Axis_LW)

        plot([-0.5 2],[0 0],'k','LineWidth',Axis_LW)

    else

        plot([0 0],[-1 25],'k','LineWidth',Axis_LW)

        plot([-1 20],[0 0],'k','LineWidth',Axis_LW)

    end

    % Trajectory
    plot(U(1:u_below),V(1:u_below),'-','LineWidth',Traj_LW,'Color',traj_col)

    plot(U(u_below:end),V(u_below:end),':','LineWidth',Traj_LW,'Color',traj_col)

    % Boundary condition
    BC_array = linspace(0,u_intercept,100);

    plot(BC_array,w_func(BC_array),'-','LineWidth',BC_LW,'Color',BC_col)

    % Intersections
    plot(u_intercept,w_func(u_intercept),'o','MarkerSize',M_size,'MarkerFaceColor',Our_BC_col, ...
        'MarkerEdgeColor',Our_BC_col)

    if c > 0

        plot(0,v_val_FS,'o','MarkerSize',M_size,'MarkerFaceColor',FS_BC_col, ...
            'MarkerEdgeColor',FS_BC_col)

        plot(0.5,v_val_Maud,'o','MarkerSize',M_size,'MarkerFaceColor',Maud_BC_col, ...
            'MarkerEdgeColor',Maud_BC_col)

    end

    % Equilibria
    plot(1,0,'ko','MarkerSize',M_size,'MarkerFaceColor','k')

    plot(0,0,'ko','MarkerSize',M_size,'MarkerFaceColor','k')

    % Axes Limits and Labels
    if c > 0

        xlim([-0.5 1.2])
        ylim([-0.5 0.5])

    else

        xlim([-1 8])
        ylim([-1 25])

    end

    xlabel('$U(z)$','Interpreter','latex')
    ylabel('$W(z)$','Interpreter','latex')

    set(gca,'FontSize',Axis_FS,'TickLabelInterpreter','latex')

    box on
    grid off

end

%%  Density profile plotting function
function plot_profile(c,tile_num)

    % Plot Styles
    Traj_LW  = evalin('base','Traj_LW');
    Axis_FS  = evalin('base','Axis_FS');
    Axis_LW  = evalin('base','Axis_LW');
    M_size   = evalin('base','M_size');

    traj_col     = evalin('base','traj_col');
    Our_BC_col   = evalin('base','Our_BC_col');
    FS_BC_col    = evalin('base','FS_BC_col');
    Maud_BC_col  = evalin('base','Maud_BC_col');

    % Solve the ODE system
    sol = solve_phase_plane(c);

    t = sol.x';

    U = sol.y(1,:)';
    V = sol.y(2,:)';

    % Phase-plane trajectory
    w_func = @(u) -c*u;

    [u_intercept,~] = ...
        Find_ustar_func(c,U,V,w_func,[]);

    u_shift = interp1(U,t,u_intercept);
    indx = find((t-u_shift) > 0,1);

    % Tile
    nexttile(tile_num)
    hold on

    % Axis
    plot([-20 20],[0 0],'k','LineWidth',Axis_LW)

    % Profile
    plot(t(1:indx)-u_shift,U(1:indx),'Color',traj_col,'LineWidth',Traj_LW)

   plot(t(indx:end)-u_shift,U(indx:end),':','Color',traj_col,'LineWidth',Traj_LW)

    fplot(@(z) 1,[-100 t(1)-u_shift],'Color',traj_col,'LineWidth',Traj_LW)

    % Main intercept
    plot(0,u_intercept,'o','MarkerSize',M_size,'MarkerFaceColor',Our_BC_col, ...
        'MarkerEdgeColor',Our_BC_col)

    % Comparison to Fisher-Stefan and modified Fisher-Stefan models
    if c > 0

        FS_z = interp1(U(1:200),t(1:200),0);
        Maud_z = interp1(U,t,0.5);

        plot(FS_z-u_shift,0,'o','MarkerSize',M_size,'MarkerFaceColor',FS_BC_col, ...
            'MarkerEdgeColor',FS_BC_col)

        plot(Maud_z-u_shift,0.5,'o','MarkerSize',M_size,'MarkerFaceColor',Maud_BC_col, ...
            'MarkerEdgeColor',Maud_BC_col)

        xlim([-10 10])
        ylim([-0.5 1.2])

    else
        xlim([-5 2])
        ylim([0 8])
    end

    % Labels
    xlabel('$z$','Interpreter','latex')
    ylabel('$U(z)$','Interpreter','latex')

    set(gca,'FontSize',Axis_FS,'TickLabelInterpreter','latex')

    box on
    grid off

end