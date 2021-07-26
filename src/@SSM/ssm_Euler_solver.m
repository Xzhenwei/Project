function p=ssm_Euler_solver(obj, N,T0,PSD,f,m,Wnode,R0) %% change the name 
    % l is the dim of manifold, n is dim of system, f= dim of filter;
maxiter=1000; tol=1e-9;

detT=T0/N;
n=obj.System.n;
M=PSD.Mz;  C=PSD.Cz;   K=PSD.Kz;   S=PSD.S; sigma=sqrt(S*2*pi); %variance

Gs=[PSD.G;zeros(n,f)]; %Gs=[Gs,zeros(2*n,1)];
z=zeros(f,N+1); v=zeros(f,N+1); p=zeros(m,N+1);
q=[z;v;p]; 
dW=zeros(f+f+m,1); 

display = obj.ssmSEulerTimeDisp;

for i=1:N

    detu=sigma*randn*sqrt(detT);  dW(2*f)=detu; 
    error1 = 1e8;
    iter = 0;
    
    while error1 > tol
        iter = iter+1; % update iteration
        z(:,i+1)=q(1:f,i+1);v(:,i+1)=q(f+1:2*f,i+1);p(:,i+1)=q(2*f+1:end,i+1);
    
        F=[z(:,i+1)-v(:,i+1)*detT-z(:,i);...
        (M+C*detT)*v(:,i+1)+K*z(:,i+1)*detT-M*v(:,i);...
        p(:,i+1)-expand_coefficients(R0,m, p(:,i+1))*detT-Wnode*Gs*v(:,i+1)*detT-p(:,i)]...
        -dW;
    
        Jp = compute_J_R0(R0,m, p(:,i+1));
       
        J= [eye(f),-eye(f)*detT ,zeros(f,m);K*detT,M+C*detT,zeros(f,m);zeros(m,f),-Wnode*Gs*detT,eye(m)-Jp*detT];
        y = -J\F;  % solve the linear equations
        q(:,i+1) = q(:,i+1) + y;
 
        % calculate errors
        error1 = sqrt(sum(y.^2));
    
    % break computation if maximum number of iterations is exceeded
        if iter == maxiter
            warning('Maximum number of iterations (%i) exceeded, solver may not have converged.', maxiter);
            break;
        end 
        
    end

    
    %% update

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

                J =  J + tenmat(T{j},max(size(T{j}))).data;
            else

                J =  J + expand_tensor_derivative(T{j},p);
            end
    end

end

