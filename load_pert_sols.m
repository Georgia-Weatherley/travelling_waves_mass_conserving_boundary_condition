function [W_pert_neg_c,W_pert_small_c,c_equ_neg_c,c_equ_small_c] = load_pert_sols
% function to load the analytic expressions for all perturbation solutions 

% Expansion for W(U) - very neg c 
W_0 = @(U) 1 - U; 
W_1 = @(U) 0; 
W_2 = @(U) -(U.^2)/2 + 1/2; 
W_3 = @(U) 0;
W_4 = @(U) (U.^2)/4 + (U.^3)/6 - 10/24; 

W_pert_neg_c = @(U,c) c*W_0(U) + W_1(U) + (1/c)*W_2(U) + ...
    (1/c^2)*W_3(U) + (1/c^3)*W_4(U);

% Expansion for W(U) - small c 
W_0 = @(U) -(1-U).*sqrt(((2*U)+1)/3);
W_1 = @(U) ((-(U-2).*(1+2*U).^(3/2)) - sqrt(27))./(5*(U-1).*sqrt(1+2.*U));

W_pert_small_c = @(U,c) W_0(U) + c*W_1(U);

% Expansion for u* vs c
c_equ_small_c = @(U) -(U-1) + (1/6)*(U-1).^2;
c_equ_neg_c = @(U) -U/sqrt(2) + 1/(3*sqrt(2)); 

