function linear_analytic=compute_linear_PSD(obj,omega,PSD)
% linear_analytic. This function computes the correspoding linear response
% PSD with the input PSD along with its frequency domain vector omega.

n = obj.n;
M = obj.M; C = obj.C; K = obj.K;
forcingdof = obj.forcingdof;

Phi_F = zeros(n,n);
linear_analytic = zeros(n,length(omega));

for j=1:length(omega)
 Hw = inv(-omega(j)^2*M+1i*omega(j)*C+K);
 Phi_F(forcingdof,forcingdof) = PSD(j);
 L = (-omega(j)^2*M+1i*omega(j)*C+K)\Phi_F * Hw.';
 
 linear_analytic(:,j) = abs(diag(L));
end