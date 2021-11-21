function p=ssm_Implicit_solver(obj, N,T0,PSD,f,m,Wnode,R0,inputForcingType) %% change the name 
    % l is the dim of manifold, n is dim of system, assmue 1-dim filter;

maxiter=1000; tol=1e-10;

detT=T0/N;
n=obj.System.n;
M=PSD.Mz;  C=PSD.Cz;   K=PSD.Kz;   S=PSD.S; sigma=sqrt(S*2*pi); %variance

Gs=[PSD.G;zeros(n,f)];
z=zeros(f,N+1); v=zeros(f,N+1); p=zeros(m,N+1);

dW=zeros(f,1); 

display = obj.ssmSEulerTimeDisp;

for i=1:N

    detu=sigma*randn*sqrt(detT);  dW(end)=detu; 
    z(:,i+1)=z(:,i)+detT*v(:,i);
    v(:,i+1)=v(:,i)-M\(detT*(C*v(:,i)+K*z(:,i))+dW);
    
    error1 = 1e8;
    iter = 0;
    switch inputForcingType
        case 'disp'
            Zf = z(:,i+1); % taking displacement
        case 'vel'
            Zf = v(:,i+1); % taking velocity
    end
    
    while error1 > tol
        F = p(:,i+1)-p(:,i)-expand_coefficients(R0,m, p(:,i+1))*detT-Wnode*Gs*Zf*detT;
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
    if norm(z(:,i+1))>1e10
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
