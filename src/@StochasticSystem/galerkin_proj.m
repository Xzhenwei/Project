function [w,Gz] = galerkin_proj(obj, V, PSD, SDEmethod, PSDpair)
assert(size(V,2)==1,'The projection vector is not 1-dim!')
Mt = V'*obj.M*V; Ct = V'*obj.C*V; Kt = V'*obj.K*V;
tol = 1e-8; maxiter = 1000;
inputForcingType = obj.inputForcingType;
dampingRatio = Ct/Mt/2/sqrt(Kt/Mt);
disp(['The damping ratio is: ', num2str(dampingRatio)])
% [~, e, ~] = eigs(sparse(obj.K),sparse(obj.M),2,'smallestabs');
% e = sqrt(e(1,1)); lambda = sqrt(Kt/Mt);
% disp(num2str(e-lambda))
switch SDEmethod
    case "filter ImplicitMidPoint"
        N = obj.nPoints;
        T0 = obj.timeSpan;
        detT = T0/N; t_span = 0:detT:T0;

        Mz=PSD.Mz;  Cz=PSD.Cz;   Kz=PSD.Kz;   S=PSD.S; G = PSD.G; sigma = sqrt(S*2*pi); %variance
        m = size(Mz,1);

        M = [Mz, sparse(m,1); sparse(1,m), Mt];
        switch inputForcingType
            case 'disp'
                C = [Cz, sparse(m,1); sparse(1,m), Ct]; % filter with nonzero mean
                K = [Kz, sparse(m,1); -V'*G, Kt];
            case 'vel'
                C = [Cz, zeros(m,1); -V'*G, Ct]; % filter with zero mean
                K = [Kz, zeros(m,1); zeros(1,m), Kt];
        end
        q = zeros(m+1,N+1); qd = zeros(m+1,N+1); 
        if obj.galerkinHeun
            for i = 1:N
                detu = sigma*randn*sqrt(detT); dW = zeros(m+1,1); dW(m) = detu;
                fnlxhat = ~obj.linear*[sparse(m,1); V'*obj.compute_fnl(V*q(m+1:end,i),V*qd(m+1:end,i))];
                qhat = q(:,i)+detT*qd(:,i);
                qdhat = qd(:,i)-M\(detT*(C*qd(:,i)+K*q(:,i)+fnlxhat)) + M\dW;
                qbar = (q(:,i) + qhat)/2;
                qdbar = (qd(:,i) + qdhat)/2;
                
                fnlx = ~obj.linear*[sparse(m,1); V'*obj.compute_fnl(V*qbar(m+1:end),V*qdbar(m+1:end))];
                q(:,i+1) = q(:,i)+detT*qdbar;
                qd(:,i+1) = qd(:,i)-M\(detT*(C*qdbar+K*qbar+fnlx)) + M\dW;
                if norm(q(:,i+1))>1e20
                    error('narrowing time step!')
                end
                if obj.sdeImpTimeDisp
                    disp(['time integration completed: ', num2str(i/N*100), '%']) 
                end
                
            end
        else

            for i = 1:N
            %% predict initial iteration point
            

            error1 = 1e8;
            qhat = q(:,i); qdhat = qd(:,i);
            iter = 0;
        %     tStart = tic;

            while error1 > tol
                iter = iter+1; % update iteration

                fnlx = ~obj.linear*[sparse(m,1); V'*obj.compute_fnl(V*qhat(m+1:end),V*qdhat(m+1:end))];
                f = (M+1/2*C*detT)*(qhat-q(:,i))+1/4*K*qhat*detT^2+1/4*fnlx*detT^2-1/2*detT*(M*qd(:,i)+1/2*dW);
                %% evaluate jacobian
                dfnldq = ~obj.linear*blkdiag(zeros(m,m),V' *obj.compute_dfnldx(V*qhat(m+1:end),V*qdhat(m+1:end))* V );

                J = M+1/2*C*detT+1/4*detT^2*(K+dfnldq);

                y = -J\f;  % solve the linear equations
                qhat = qhat + y(1:m+1);

                % calculate errors
                error1 = sqrt(sum(y.^2));

            % break computation if maximum number of iterations is exceeded
                if iter == maxiter
                    warning('Maximum number of iterations (%i) exceeded, solver may not have converged.', maxiter);
                    break;
                end

            end
            fnl_hat = ~obj.linear*[sparse(m,1); V'*obj.compute_fnl(V*qhat(m+1:end),V*qdhat(m+1:end))];
            qdhat = (M+1/2*C*detT)\(M*qd(:,i)-1/2*(K*qhat+fnl_hat)*detT+1/2*dW);

        %     tEnd = toc(tStart)

            %% update
            q(:,i+1) = q(:,i)+detT*qdhat;
            xhat = qhat(m+1:end);xdhat=qdhat(m+1:end);
            fnl_hat = ~obj.linear*[sparse(m,1); V'*obj.compute_fnl(V*xhat,V*xdhat)];
            qd(:,i+1) = qd(:,i)-M\(C*qdhat+K*qhat+fnl_hat)*detT+M\dW;
                if obj.sdeImpTimeDisp
                disp(['time integration completed: ', num2str(i/N*100), '%'])   
                end
            end
        end
        
        X = q(m+1:end,:); x = V*X ;
        
    case "Newmark"
        %%% when specifying a given 
        obj.Fsto = obj.generate_stochastic();
        tol = 1e-8;
        N = obj.nPoints;
        T0 = obj.timeSpan;
        detT = T0/N; t_span = 0:detT:T0;

        M =  Mt; C= Ct; K=Kt;
        q = zeros(1,N+1); qd = zeros(1,N+1); 
        %%
        w = obj.samplePSD(2,:);
        nOmega = length(w); linear_galerkin = zeros(1, nOmega);
        for j = 1:nOmega 
                PSDmatrix = zeros(obj.n,obj.n) ;
                PSDmatrix(end,end) = obj.samplePSD(1,j);
                Zj = V'* PSDmatrix *V;
%                 Zj = V(2)^2 * obj.samplePSD(1,j) ;
                Hw = 1/(-w(j)^2*M+1i*w(j)*C+K);
                Z_full = abs(Hw)^2*Zj; 

                linear_galerkin(j) = Z_full;
        end
        numerator = 1;
        denominator = [M,C,K];
        sys = tf(numerator,denominator);
        bode(sys)
        figure
        plot(w,20*log(linear_galerkin))
        xlim([2 25])

        %%
        %Heun on original system

        for i = 2:N+1            
            t = (i-1)*detT;
            qhat = q(i-1) + detT*qd(i-1);
            Fext = V'*obj.compute_fstochastic(t);
            Fqxn = V'*obj.compute_fnl(V*q(i-1),V*qd(i-1));
            qdhat = qd(i-1) - M\C*qd(i-1)*detT - M\K*q(i-1)*detT-M\Fqxn*detT + M\Fext*detT;
            
            Fqhat = V'*obj.compute_fnl(V*qhat,V*qdhat);
            q(i) = q(i-1) + detT*(qd(i-1)+qdhat)/2;
            qd(i) = qd(i-1) - M\C*(qd(i-1)+qdhat)/2*detT-...
                M\(Fqxn+Fqhat)/2*detT - M\K*(q(i-1)+qhat)/2*detT + M\Fext*detT;
                if norm(q(:,i)) > 1e10
                    error('Forward step does not converge, please decrease the time step')
                end
%                 disp(['time integration completed: ', num2str(i/N*100), '%'])            
        end

         x = V*q ;
        %     %% ODE45 solver
%                 f_p = @(t,q) Galerkin(obj,V,Mt,Ct,Kt,t,q);
%                 opts = odeset('RelTol',tol,'AbsTol',tol);
%                 [t_span,q_45] = ode45(f_p,t_span ,zeros(2,1), opts);
%                 q_45 = q_45';
%                 x = V.*q_45(1,:);
        

end
%%
                [w,Gz] = crossPSDestimator(x(PSDpair(1),:),x(PSDpair(2),:),t_span);
                

                
end

function dq = Galerkin(obj,V,Mt,Ct,Kt,t,q)
            dq1 = q(2);
            dq2 = (V'*obj.compute_fstochastic(t) - V'*obj.compute_fnl(V*q(1),V*q(2)) - Kt*q(1) - Ct*q(2))/Mt ;
            dq = [dq1;dq2];
end
function [w,Gxy]=crossPSDestimator(x,y,t)
%%% Reference: Preunont Andre, Random vibration and Spectral Analysis Chp 12.9
%%%% numebr of sampling points 
N = length(x); 
%%%% time period of generation
T0 = max(t);
%%%% sampling frequency, time spacing
f0 = 1/T0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cx = fft(x)/N; Cy = fft(y)/N; Cys = conj(Cy);

Gxy = T0*(Cx.*Cys)/2/pi;
w = (1:length(Cx))*f0*2*pi;

end