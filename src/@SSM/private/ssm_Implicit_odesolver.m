function p=ssm_Implicit_odesolver(obj, N,T0,m,Wnode,R0) %% change the name 
    % l is the dim of manifold, n is dim of system, assmue 1-dim filter;

maxiter=1000; tol=1e-10;
obj.System.Fsto = obj.System.generate_stochastic();
detT=T0/N;
n=obj.System.n;

p=zeros(m,N+1);

display = obj.ssmSEulerTimeDisp;

for i=1:N

    
    error1 = 1e8;
    iter = 0;
    Fext = [obj.System.compute_fstochastic(detT*(i-1)); sparse(n,1)];
    while error1 > tol
        
        F = p(:,i+1)-p(:,i)-expand_coefficients(R0,m, p(:,i+1))*detT-Wnode*Fext*detT;
        Jp = compute_J_R0(R0,m, p(:,i+1))*detT;
        J= eye(m) - Jp;
        y = -J\F;  % solve the linear equations
        p(:,i+1) = p(:,i+1) + y;
 
        % calculate errors
        error1 = sqrt(sum(y.^2));
    
    % break computation if maximum number of iterations is exceeded
        if iter == maxiter
            warning('Maximum number of iterations (%i) exceeded, solver may not have converged.', maxiter);
            break;
        end 
    end
    %% update
    if norm(p(:,i+1))>1e10
        error('narrowing time step')
    end
    if display
        disp(['time integration completed: ', num2str(i/N*100), '%']) 
    end
end

end
    

function S = expand_coefficients(S0,n, p)
% m = numel(p);
S = zeros(n,1);

% expand autonomous coefficients
    for k = 1:length(S0)
        S =  S + expand_multiindex(S0{k},p);
    end
    
end

function J = compute_J_R0(R0,n, p)
% COMPUTE_J_R0 This function computes the Jacobian of  nonlinear R0
% with respect to the reduced coordinate p in a first-order
% 
J = zeros(n,n); T={};l=0;

% transfer multi-index to tensor
    for k = 1:length(R0)
        if ~isempty(R0{k}.coeffs)
            T{l+1} = multi_index_to_tensor(R0{k}.coeffs, R0{k}.ind);
            l=l+1;
        end
    end
    for j=1:l
            if length(size(T{j}))==length(p)
                Tensor2Mat = tenmat(T{j},max(size(T{j})));
                J =  J + Tensor2Mat.data;
            else

                J =  J + expand_tensor_derivative(T{j},p);
            end
    end

end
