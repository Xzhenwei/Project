function [M,C,K,fnl] = build_model(n,m,c,k,kappa2,kappa3)

[K,C,~,~] = assemble_global_coefficients(k,kappa2,kappa3,c,n);
%%
% Dirichlet boundary conditions

K = K(2:n+1,2:n+1); K(end,end)=K(end,end)/2;
C = C(2:n+1,2:n+1); C(end,end)=C(end,end)/2;
M = m*speye(n,n);

f2 = sptensor([n,n,n]); % coeff of quadratic nonlinearity/ 0 
f3 = sptensor([n,n,n,n]); % coeff of cubic nonlinearity / single cubic spring
%% 
% Adding cubic spring to first node of the beam 
f3(1,1,1,1) = kappa3;

fnl = {f2, f3};