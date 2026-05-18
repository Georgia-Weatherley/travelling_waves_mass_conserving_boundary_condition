function sol = solve_phase_plane(c)
% This function solves the ODE system for a given wave speed c 

    % Linearisation about (1,0)
    A = [0 1; 1 -c];
    [V,D] = eig(A);
    eigvals = diag(D);

    %% Choose unstable/stable eigenvector
    if c >= 0
        idx = eigvals > 0;
    else
        idx = eigvals < 0;
    end

    vec = V(:,idx);

    % Small perturbation from equilibrium
    eps = 1e-4;
    y0 = [1;0] + eps*vec;

    % Special case c = 0
    if c == 0
        y0 = [0.99; -5e-4];
    end

    %%
    % Integration interval
    zspan = [0 200];

    % Solver options
    opts = odeset( ...
        'RelTol',1e-9, ...
        'AbsTol',1e-11);

    % Solve
    sol = ode45(@(z,y) phase_ode(z,y,c), ...
                zspan, y0, opts);

end

function dydz = phase_ode(~,y,c)
% function that defines the ODE system 

    U = y(1); V = y(2);
    dydz = [V; -c*V - U*(1-U)];

end