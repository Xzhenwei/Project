function [M,C,K,fnl,k1,c1] = build_model(kappa, nElements)

[M,C,K,k1,c1]=L_QuarterCar_Model(nElements);

n = length(M);
% forcingDof = n;



%% can be transfered to QC model
% first-order tensors

f2 = sptensor([n,n,n]); % coeff of quadratic nonlinearity/ 0 
f3 = sptensor([n,n,n,n]); % coeff of cubic nonlinearity / single cubic spring
%% 
% Adding cubic spring to end node of the cart
dof =  n;
for j = dof
    f3(j,j,j,j) = kappa;
end
fnl = {f2,f3};
