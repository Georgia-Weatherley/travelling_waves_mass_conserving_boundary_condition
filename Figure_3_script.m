clear; clc; close all;
% This script reproduces Figure 3 from the manuscript. 
% It calls the functions: 
% solve_phase_plane.m 
% Find_ustar_func.m 

% It also generates and saves the wave speed curve (u* vs c) as
% "wave_speed_curve.mat"

%%  Figure setup and plot styles
Traj_LW = 6;
M_size  = 19;
Axis_LW = 1.5; 
Axis_FS = 30;
traj_col = [0.2431    0.2980    0.3686];

map = 1/255*[77,0,75;
             129,15,124;
             136,65,157;
             140,107,177;
             140,150,198;
             158,188,218;
             191,211,230;
             166,189,219;
             103,169,207;
             54,144,192;
             2,129,138;
             1,108,89;
             1,70,54];

figure
tiledlayout(2,2)

%%  Positive c values
nexttile(1), hold on
c_array = fliplr([0.25 0.5 1 1.5 2 10]);

alpha_array = zeros(size(c_array));
u_intercept_array = zeros(size(c_array));

C = flipud(1/255*[140,150,198;
                  140,107,177;
                  136,65,157;
                  129,15,124;
                  82,0,78;
                  41,0,39]);

plot([0 0],[-0.5 2],'k','LineWidth',Axis_LW)
plot([-0.5 2],[0 0],'k','LineWidth',Axis_LW)

posc_u_intercept_array = zeros(size(c_array));

for i = 1:length(c_array)

    c = c_array(i);
    w_func = @(u) -c*u;

    sol = solve_phase_plane(c);

    U = sol.y(1,:)';
    V = sol.y(2,:)';

    if c > 0 && c < 2

        [u_intercept, u_below] = Find_ustar_func(c,U,V,w_func,[]);

        alpha_array(i) = c / u_intercept;
        u_intercept_array(i) = u_intercept;
        posc_u_intercept_array(i) = u_intercept;

        txt = ['$c = $' num2str(c)];

        plot(U(1:u_below),V(1:u_below),'LineWidth',Traj_LW,'Color',C(i,:), ...
            'DisplayName',txt)

        plot(U(u_below:end),V(u_below:end),':','LineWidth',2,'Color',[C(i,:) 0.5], ...
            'HandleVisibility','off')

        plot(u_intercept,w_func(u_intercept),'o','MarkerSize',M_size,'MarkerFaceColor',C(i,:), ...
            'MarkerEdgeColor',C(i,:),'HandleVisibility','off')

    else

        txt = ['$c = $' num2str(c)];

        plot(U,V,'-','LineWidth',Traj_LW,'Color',C(i,:),'DisplayName',txt)
    end
end

plot(1,0,'ko','MarkerFaceColor','k','HandleVisibility','off')
plot(0,0,'ko','MarkerFaceColor','k','HandleVisibility','off')

xlim([-0.4 1.1])
ylim([-0.5 0.4])

xlabel('$U(z)$','Interpreter','latex')
ylabel('$W(z)$','Interpreter','latex')

set(gca,'FontSize',Axis_FS,'TickLabelInterpreter','latex')
box on; grid off

pos_c_array = c_array;
pos_col = C;

%%  Negative c values

nexttile(2), hold on

c_array = fliplr([-0.5 -1 -2 -3 -4 -5]);

u_intercept_array = zeros(size(c_array));

C = flipud(1/255*[166,189,219;
                  103,169,207;
                  54,144,192;
                  2,129,138;
                  1,108,89;
                  1,70,54]);

plot([0 0],[-5 50],'k','LineWidth',Axis_LW)
plot([-5 50],[0 0],'k','LineWidth',Axis_LW)

for i = 1:length(c_array)

    c = c_array(i);
    w_func = @(u) -c*u;

    sol = solve_phase_plane(c);

    U = sol.y(1,:)';
    V = sol.y(2,:)';

    [u_intercept, u_below] = Find_ustar_func(c,U,V,w_func,[]);

    u_intercept_array(i) = u_intercept;

    txt = ['$c = $' num2str(c)];

    plot(U(1:u_below),V(1:u_below),'-','LineWidth',Traj_LW,'Color',C(i,:), ...
        'DisplayName',txt)

    plot(U(u_below:end),V(u_below:end),':','LineWidth',2,'Color',[C(i,:) 0.5], ...
        'HandleVisibility','off')

    plot(u_intercept,w_func(u_intercept),'o','MarkerSize',M_size,'MarkerFaceColor',C(i,:), ...
        'MarkerEdgeColor',C(i,:),'HandleVisibility','off')

end

plot(1,0,'ko','MarkerFaceColor','k','HandleVisibility','off')
plot(0,0,'ko','MarkerFaceColor','k','HandleVisibility','off')

xlim([0 8])
ylim([0 40])

xlabel('$U(z)$','Interpreter','latex')
ylabel('$W(z)$','Interpreter','latex')

set(gca,'FontSize',Axis_FS,'TickLabelInterpreter','latex')
box on; grid off

negc_u_intercept_array = u_intercept_array;
neg_c_array = c_array;
neg_col = C;

%%  Wave speed (c vs u*) curve

figure, hold on

c_array = linspace(-10,1.98,200);
u_star_PP = zeros(size(c_array));

for i = 1:length(c_array)

    c = c_array(i);
    w_func = @(u) -c*u;

    sol = solve_phase_plane(c);

    U = sol.y(1,:)';
    V = sol.y(2,:)';

    u_star_PP(i) = Find_ustar_func(c,U,V,w_func,[]);
end

plot(u_star_PP,c_array,'-','Color',traj_col,'LineWidth',Traj_LW)

% Positive markers
for ii = 1:length(pos_c_array)
    if pos_c_array(ii) < 2
        plot(posc_u_intercept_array(ii),pos_c_array(ii),'o', ...
            'MarkerSize',M_size, ...
            'MarkerFaceColor',pos_col(ii,:), ...
            'MarkerEdgeColor',pos_col(ii,:), ...
            'DisplayName',num2str(pos_c_array(ii)))
    end
end

% Negative markers
neg_c_array = fliplr(neg_c_array);
negc_u_intercept_array = fliplr(negc_u_intercept_array);
neg_col = flipud(neg_col);

for jj = 1:length(neg_c_array)
    plot(negc_u_intercept_array(jj),neg_c_array(jj),'o', ...
        'MarkerSize',M_size, ...
        'MarkerFaceColor',neg_col(jj,:), ...
        'MarkerEdgeColor',neg_col(jj,:), ...
        'DisplayName',num2str(neg_c_array(jj)))
end

% Plot style
xlim([0 8])
ylim([-6 2])

legend('NumColumns',10,'Location','southoutside','Interpreter','latex')
xlabel('$U^{*}$','Interpreter','latex')
ylabel('$c$','Interpreter','latex')

set(gca,'FontSize',Axis_FS,'TickLabelInterpreter','latex')
box on

save("wave_speed_curve.mat", "u_star_PP", "c_array")