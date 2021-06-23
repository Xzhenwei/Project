function [M,C,K,fnl,f_0] = build_model(kappa, nElements)

[M,C,K]=L_QuarterCar_Model(nElements);

n = length(M);
forcing_dof = n;

f_0 = sparse(n,1);
f_0(forcing_dof) = 1;


%% can be transfered to QC model
% first-order tensors

f2 = sptensor([n,n,n]); % coeff of quadratic nonlinearity/ 0 
f3 = sptensor([n,n,n,n]); % coeff of cubic nonlinearity / single cubic spring
%% 
% Adding cubic spring to end node of the beam 
dof =  n;
for j = dof
    f3(j,j,j,j) = kappa;
end
fnl = {f2,f3};
