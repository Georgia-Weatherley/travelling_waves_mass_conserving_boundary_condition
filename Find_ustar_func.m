function [u_star,u_below] = Find_ustar_func(c,U,V,w_func,uval)
% Identifies where the boundary function intersects the solution
% trajectory to determine the front density u*

if isempty(uval) 
    if c > 0
        mask = find(U < 0);
        New_V = V(1:mask(1));
        New_U = U(1:mask(1));
        shifted_V = New_V - w_func(New_U);
        mask2 = find(shifted_V < 0);
        u_below = mask2(1); u_above = mask2(1)-1;
        p = polyfit(New_U([u_below,u_above]), shifted_V([u_below,u_above]), 1);
        u_star = -p(2)/p(1);
       
    elseif c == 0
        u_star = 1;
        u_below = [];

    else % c < 0
        shifted_V = w_func(U) - V;
        mask2 = find(shifted_V < 0);
        u_below = mask2(1); u_above = mask2(1)-1;
        p = polyfit(U([u_below,u_above]), shifted_V([u_below,u_above]), 1);
        u_star = -p(2)/p(1);

    end

else
       uval_array = uval*ones(length(U),1);
        diff = U - uval_array;
        mask = find(diff < 0);
        u_below = mask(1);
        u_above = u_below - 1;
        p = polyfit(U([u_below, u_above]), V([u_below, u_above]), 1);
        u_star = polyval(p, uval);
end