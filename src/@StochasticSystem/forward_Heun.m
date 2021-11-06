function [X,V] = forward_Heun(obj,N,T0,PSD)
    n = obj.n;
    Mz = PSD.Mz;
    Cz = PSD.Cz;
    Kz = PSD.Kz;
    S = PSD.S;
    G = PSD.G;   G(obj.forcingdof,1) = 1;
    m = length(Mz);       
    detT = T0/N;  
    
    sigma = sqrt(S*2*pi); %variance
    M = [Mz, sparse(m,obj.n); sparse(obj.n,m), obj.M];
    C = [Cz, sparse(m,obj.n); sparse(obj.n,m), obj.C];
    K = [Kz, sparse(m,obj.n); -G, obj.K];
    q = zeros(m+n,N+1); qd=zeros(m+n,N+1);

        isnonlinear = ~obj.linear;
        %Heun on original system

        for i = 2:N+1            
            detu = sigma*randn*sqrt(detT);
            detF = [zeros(m-1,1); detu]; dF = [detF;zeros(n,1)];

            qhat = q(:,i-1) + detT*qd(:,i-1);
            Fqxn = q(:,i-1) + isnonlinear*[zeros(m,1);compute_fnl(obj,q(m+1:end,i-1),qd(m+1:end,i-1))];
            qdhat = qd(:,i-1) - M\C*qd(:,i-1)*detT - M\K*detT*Fqxn + M\dF;
            
            Fqhat = qhat + isnonlinear * [zeros(m,1);compute_fnl(obj,qhat(m+1:end),qdhat(m+1:end))];
            q(:,i) = q(:,i-1) + detT*(qd(:,i-1)+qdhat)/2;
            qd(:,i) = qd(:,i-1) - M\C*(qd(:,i-1)+qdhat)/2*detT-...
                M\K*(Fqxn+Fqhat)/2*detT + M\dF;
                if isnan(q(:,i))
                    error('Forward step does not converge, please decrease the time step')
                end
%                 disp(['time integration completed: ', num2str(i/N*100), '%'])            
        end
X = q(m+1:end,:); V = qd(m+1:end,:); 
end