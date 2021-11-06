function [X,V]=implicit_Mid_Point(obj,N,T0,PSD)
    maxiter = 10000; tol = 1e-12;

    n = obj.n;
    Mz = PSD.Mz;
    Cz = PSD.Cz;
    Kz = PSD.Kz;
    S = PSD.S;
    G = PSD.G;   
    m = length(Mz);
    detT = T0/N;  
    
    sigma = sqrt(S*2*pi); %variance
    M = [Mz, sparse(m,n); sparse(n,m), obj.M];
%     C = [Cz, sparse(m,n); sparse(n,m), obj.C]; % filter with nonzero mean
%     K = [Kz, sparse(m,n); -G, obj.K];

    C = [Cz, sparse(m,n); -G, obj.C]; % filter with zero mean
    K = [Kz, sparse(m,n); sparse(n,m), obj.K];
    
    q = zeros(m+n,N+1); qd = zeros(m+n,N+1); 

%% below only for fxnl(x)
    for i = 1:N
    %% predict initial iteration point
    detu = sigma*randn*sqrt(detT); dW = zeros(m+n,1); dW(m) = detu;
    
    error1 = 1e8;
    qhat = q(:,i); qdhat = qd(:,i);
    iter = 0;
%     tStart = tic;

    while error1 > tol
        iter = iter+1; % update iteration
        
        fnlx = ~obj.linear*[sparse(m,1); obj.compute_fnl(qhat(m+1:end),qdhat(m+1:end))];
        f = (M+1/2*C*detT)*(qhat-q(:,i))+1/4*K*qhat*detT^2+1/4*fnlx*detT^2-1/2*detT*(M*qd(:,i)+1/2*dW);
        %% evaluate jacobian
        dfnldq = ~obj.linear*blkdiag(zeros(m,m),obj.compute_dfnldx(qhat(m+1:end),qdhat(m+1:end)));

        J = M+1/2*C*detT+1/4*detT^2*(K+dfnldq);

        y = -J\f;  % solve the linear equations
        qhat = qhat + y(1:m+n);

        % calculate errors
        error1 = sqrt(sum(y.^2));
    
    % break computation if maximum number of iterations is exceeded
        if iter == maxiter
            warning('Maximum number of iterations (%i) exceeded, solver may not have converged.', maxiter);
            break;
        end
        
    end
    fnl_hat = ~obj.linear*[sparse(m,1); obj.compute_fnl(qhat(m+1:end),qdhat(m+1:end))];
    qdhat = (M+1/2*C*detT)\(M*qd(:,i)-1/2*(K*qhat+fnl_hat)*detT+1/2*dW);

%     tEnd = toc(tStart)
    
    %% update
    q(:,i+1) = q(:,i)+detT*qdhat;
    xhat = qhat(m+1:end);xdhat=qdhat(m+1:end);
    fnl_hat = ~obj.linear*[sparse(m,1); obj.compute_fnl(xhat,xdhat)];
    qd(:,i+1) = qd(:,i)-M\(C*qdhat+K*qhat+fnl_hat)*detT+M\dW;
        if obj.sdeImpTimeDisp
        disp(['time integration completed: ', num2str(i/N*100), '%'])   
        end
    end
    
X = q(m+1:end,:); V = qd(m+1:end,:);
end
