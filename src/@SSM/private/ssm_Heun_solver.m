function p=ssm_Heun_solver(obj, N,T0,PSD,f,m,Wnode,R0, inputForcingType) %% change the name 
    % l is the dim of manifold, n is dim of system, f= dim of filter;

detT=T0/N;
n=obj.System.n;
M=PSD.Mz;  C=PSD.Cz;   K=PSD.Kz;   S=PSD.S; sigma=sqrt(S*2*pi); %variance

Gs=[PSD.G;zeros(n,f)]; %Gs=[Gs,zeros(2*n,1)];
z=zeros(f,N+1); v=zeros(f,N+1); p=zeros(m,N+1);

dW=zeros(f+f+m,1); 

display = obj.ssmSEulerTimeDisp;

for i=1:N

    detu = sigma*randn*sqrt(detT);  dW(2*f)=detu; 
    zhat = z(:,i)+detT*v(:,i);
    vhat = v(:,i)-M\(detT*(C*v(:,i)+K*z(:,i))+detu);
    phat = p(:,i)+expand_coefficients(R0,m, p(:,i))*detT+Wnode*Gs*v(:,i)*detT;
    
    zbar = (z(:,i) + zhat)/2;
    vbar = (v(:,i) + vhat)/2;
    pbar = (p(:,i) + phat)/2;
    
    z(:,i+1)=z(:,i)+detT*vbar;
    v(:,i+1)=v(:,i)-M\(detT*(C*vbar+K*zbar)+detu); 
    switch inputForcingType
        case 'disp'
            p(:,i+1)=p(:,i)+expand_coefficients(R0,m, pbar)*detT+Wnode*Gs*zbar*detT; % taking disp
        case 'vel'
            p(:,i+1)=p(:,i)+expand_coefficients(R0,m, pbar)*detT+Wnode*Gs*vbar*detT; % taking vel
    end
    
    %% update
    if norm(p(:,1))>1e20
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


