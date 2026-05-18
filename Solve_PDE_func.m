function [t,sol,rho] = Solve_PDE_func(a,b,L,N,init_val,t_start,t_stop,t_spacing,u0_array,specialIC)

% Solves the model: Fisher-KPP governing equation, mass conserving moving
% boundary with a linear boundary velocity function 

drho = 1/(N-1); 
loggrow = 1; 

u0 = zeros(N+1,1);
rho = zeros(N,1);
u0(N+1) = L; 
f_of_u = @(u_np1) (a*u_np1)+b;

for i = 1:N
    rho(i) = (i-1)*drho;
end
if specialIC == 0
    u0(1:N) = init_val;
else
    u0(1:N) = u0_array;

end

Jcbn = zeros(N+1,N+1);
diag = logical(eye(size(Jcbn)));
upper_diag = logical([diag(2:end,:); zeros(1,N+1)]);
lower_diag = logical([zeros(1,N+1); diag(1:end-1,:)]);
Jcbn(N-1:N,:) = 0;
Jcbn(diag) = 1;
Jcbn(upper_diag) = 1;
Jcbn(lower_diag) = 1;
Jcbn(:,N) = 1;
Jcbn(:,N+1) = 1; 

[row,col,val] = find(Jcbn == 1);
S = sparse(row,col,val,N+1,N+1); 

odefunc = @(t,u) PDE_discretisation(t,u,N,rho,drho,loggrow,f_of_u);

tgrid = linspace(t_start,t_stop,t_spacing);
[t,sol] = ode15s(odefunc,tgrid,u0,odeset('JPattern',S,'RelTol',1e-8));


function  du  = PDE_discretisation(t,u,N,rho,drho,loggrow,f_of_u)

du = zeros(N+1,1);
    du(N+1) = f_of_u(u(N));

    % i = 1
    du(1) = (1/(u(N+1)^2))*((u(1+1) - 2*u(1) + u(1+1))/drho^2) + loggrow*u(1)*(1-u(1));

    for i = 2:N-1
        du(i) = (1/(u(N+1)^2))*((u(i+1) - 2*u(i) + u(i-1))/drho^2) + du(N+1)*(rho(i)/u(N+1))*((u(i+1)-u(i-1))/(2*drho)) + loggrow*u(i)*(1-u(i));
    end
    
    % i = N
    du(N) = (1/(u(N+1)^2))*((u(N-1) - 2*drho*u(N)*u(N+1)*du(N+1) - 2*u(N) + u(N-1))/drho^2) + du(N+1)*(rho(i)/u(N+1))*(-u(N)*u(N+1)*du(N+1)) + loggrow*u(N)*(1-u(N));
   
end

end