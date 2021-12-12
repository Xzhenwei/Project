function [M,C,K,fnl] = build_model(n,m,c,k,kappa2,kappa3)

[K,C,f2,f3] = assemble_global_coefficients(k,kappa2,kappa3,c,n);
%%
% One end fixed Boundary condition

K = K(3:n+2,3:n+2); 
C = C(3:n+2,3:n+2); 
M = m*speye(n,n);

f2 = f2(3:n+2,3:n+2,3:n+2);
f3 = f3(3:n+2,3:n+2,3:n+2,3:n+2);
%% 
% Adding cubic spring to first node of the beam 

fnl = {f2, f3};

