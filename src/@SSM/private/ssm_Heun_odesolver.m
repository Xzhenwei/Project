function p=ssm_Heun_odesolver(obj, N,T0,m,Wnode,R0)
obj.System.Fsto = obj.System.generate_stochastic();
detT=T0/N;
n=obj.System.n;

p=zeros(m,N+1);

display = obj.ssmSEulerTimeDisp;
for i=1:N

    Fext = [obj.System.compute_fstochastic(detT*(i-1)); sparse(n,1)];
    phat = p(:,i) + expand_coefficients(R0,m, p(:,i))*detT + Wnode*Fext*detT;
    pbar = (p(:,i) + phat)/2;
    p(:,i+1)=p(:,i)+expand_coefficients(R0,m, pbar)*detT + Wnode*Fext*detT;
    
    if norm(p(:,i+1))>1e10
        error('narrowing time step, will not converge')
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


