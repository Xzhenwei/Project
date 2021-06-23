function [X,V]=implicit_Mid_Point(obj,N,T0,PSD)
    maxiter=10000; tol=1e-7;

    n=obj.System.n;
    Mz=PSD.Mz;
    Cz=PSD.Cz;
    Kz=PSD.Kz;
    S=PSD.S;
    G= PSD.G;   G(n,1)=1;
    m=length(Mz);       
    detT=T0/N;  %T=linspace(0,T0,N+1);
    
    sigma=sqrt(S*2*pi); %variance
    M = [Mz, sparse(m,n); sparse(n,m), obj.System.M];
    C = [Cz, sparse(m,n); sparse(n,m), obj.System.C];
    K = [Kz, sparse(m,n); -G, obj.System.K];
    
    q=zeros(m+n,N+1); qd=zeros(m+n,N+1); 
    
%     isnonlinear=~obj.linear;
%     for i=1:N
%     %% predict initial iteration point
%     detu=sigma*randn*sqrt(detT); dW=zeros(m+n,1); dW(m)=detu;
% %     q0=q(:,i)+detT*M\qd(:,i) ;
% %     qd0=qd(:,i)-M\(C*qd(:,i)+K*q(:,i)+fnl)*detT+M\dW ;
%     %% find midpoint via Newton method
%     
%     error1 = 1e8;
%     qhat=q(:,i); qdhat=qd(:,i);
%     iter = 0;
% %     tStart = tic;
% 
%     while error1 > tol
%         iter = iter+1; % update iteration
%         
%         fnl=[sparse(m,1); obj.compute_fnl(qhat(m+1:end),qdhat(m+1:end))];
%         f = [qhat-1/2*qdhat*detT-q(:,i);...
%             M*qdhat+1/2*(C*qdhat+K*qhat+fnl)*detT-qd(:,i)-1/2*dW]; % evaluate function at current point
%         %% evaluate jacobian
%         dfnldq=blkdiag(zeros(m,m),obj.compute_dfnldx(q(m+1:end,i),qd(m+1:end,i)));
%         dfnldqd=blkdiag(zeros(m,m),obj.compute_dfnldxd(q(m+1:end,i),qd(m+1:end,i)));
% 
%         J1 = [eye(n+m),-1/2*detT*eye(n+m)]; 
%         J2 =[1/2*(K+dfnldq)*detT,M+1/2*(C+dfnldqd)*detT];
%         J= vertcat(J1,J2);
% 
%         y = -J\f;  % solve the linear equations
%         qhat = qhat + y(1:m+n);
%         qdhat= qdhat + y(m+n+1:end);% move the solution
% 
%         % calculate errors
%         error1 = sqrt(sum(y.^2));
% %         error(iter) = sqrt(sum(f.^2));
%     
%     % break computation if maximum number of iterations is exceeded
%         if iter == maxiter
%             warning('Maximum number of iterations (%i) exceeded, solver may not have converged.', maxiter);
%             break;
%         end 
%         
%     end
% %     iter
% %     tEnd = toc(tStart)
%     
%     %% update
%     q(:,i+1)=q(:,i)+detT*qdhat;
%     xhat=qhat(m+1:end);xdhat=qdhat(m+1:end);
%     fnl_hat=[sparse(m,1); obj.compute_fnl(xhat,xdhat)];
%     qd(:,i+1)=qd(:,i)-M\(C*qdhat+K*qhat+fnl_hat)*detT+M\dW;
%     disp(['time integration completed: ', num2str(i/N*100), '%'])    
%     end

%% below only for fxnl(x)
    for i=1:N
    %% predict initial iteration point
    detu=sigma*randn*sqrt(detT); dW=zeros(m+n,1); dW(m)=detu;
    
    error1 = 1e8;
    qhat=q(:,i); qdhat=qd(:,i);
    iter = 0;
%     tStart = tic;

    while error1 > tol
        iter = iter+1; % update iteration
        
        fnl=[sparse(m,1); obj.System.compute_fnl(qhat(m+1:end),qdhat(m+1:end))];
        f= (M+1/2*C*detT)*(qhat-q(:,i))+1/4*K*qhat*detT^2+1/4*fnl*detT^2-1/2*detT*(M*qd(:,i)+1/2*dW);
        %% evaluate jacobian
        dfnldq=blkdiag(zeros(m,m),obj.System.compute_dfnldx(qhat(m+1:end),qdhat(m+1:end)));

        J= M+1/2*C*detT+1/4*detT^2*(K+dfnldq);

        y = -J\f;  % solve the linear equations
        qhat = qhat + y(1:m+n);

        % calculate errors
        error1 = sqrt(sum(y.^2));
%         error(iter) = sqrt(sum(f.^2));
    
    % break computation if maximum number of iterations is exceeded
        if iter == maxiter
            warning('Maximum number of iterations (%i) exceeded, solver may not have converged.', maxiter);
            break;
        end
        
    end
    fnl_hat=[sparse(m,1); obj.System.compute_fnl(qhat(m+1:end),qdhat(m+1:end))];
    qdhat=(M+1/2*C*detT)\(M*qd(:,i)-1/2*(K*qhat+fnl_hat)*detT+1/2*dW);
%     iter
%     tEnd = toc(tStart)
    
    %% update
    q(:,i+1)=q(:,i)+detT*qdhat;
    xhat=qhat(m+1:end);xdhat=qdhat(m+1:end);
    fnl_hat=[sparse(m,1); obj.System.compute_fnl(xhat,xdhat)];
    qd(:,i+1)=qd(:,i)-M\(C*qdhat+K*qhat+fnl_hat)*detT+M\dW;
    disp(['time integration completed: ', num2str(i/N*100), '%'])    
    end

X=q(m+1:end,:); V=qd(m+1:end,:);
end
