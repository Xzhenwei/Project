function S = expand_coefficients(obj,S0,U,epsilon,n, t, p)
% m = numel(p);
S = zeros(n,1);

% expand autonomous coefficients
for j = 1:length(S0)
%     S =  S + real(expand_multiindex(S0{j},p));
    S =  S + expand_multiindex(S0{j},p);
end

% add nonautonomous contribution
    S = S + epsilon * U * [compute_fstochastic(obj,t); zeros(obj.n,1)]; 
    % external forcing in 2nd sys is 10x1. 1st order sys is 20x1
end