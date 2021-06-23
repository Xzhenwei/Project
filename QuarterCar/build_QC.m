clear all
clc
n=2; kappa=0;
[M_f,C_f,K_f,fnl,f_0] = build_model(kappa,n);

nRealization=1;
T0=100; %% PSD frequency domain resolution is ~ 1/T0
N=2^14; %% control the accuracy of numerical differential equation
epsilon=1; %% forcing magnitude
[filterPSD,forcingdof,IC,stochastic_f] = build_stochasticF(n,N,epsilon);
DS = DynamicalSystem();
set(DS,'M',M_f,'C',C_f,'K',K_f,'fnl',fnl,'filterPSD',filterPSD,'stochastic_f',stochastic_f,'linear',false); 
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex');
DS.add_random_forcing(nRealization, T0, N,forcingdof);
A=DS.A;B=DS.B;

PSD=filterPSD;
    Mz=PSD.Mz;
    Cz=PSD.Cz;
    Kz=PSD.Kz;
    S=PSD.S;
    G=zeros(n,length(Mz));   G(n,1)=1;
    m=length(Mz);       
    detT=T0/N;  t=linspace(0,T0,N+1);
    
    sigma=sqrt(S*2*pi); %variance
    M = [Mz, sparse(m,n); sparse(n,m), M_f];
    C = [Cz, sparse(m,n); sparse(n,m), C_f];
    K = [Kz, sparse(m,n); -G, K_f];
    q=zeros(m+n,N+1); qd=zeros(m+n,N+1); q2=q; q2d=q;
    isnonlinear=false;
         for i=2:N+1            
            detu=sigma*randn*sqrt(detT);
            detF=[zeros(m-1,1);detu]; dF=[detF;zeros(n,1)];

            qhat=q(:,i-1)+detT*qd(:,i-1);
            Fqxn=q(:,i-1)+isnonlinear*[zeros(m,1);compute_fnl(DS,q(m+1:end,i-1),qd(m+1:end,i-1))];
            qdhat=qd(:,i-1)-M\C*qd(:,i-1)*detT-M\K*detT*Fqxn+M\dF;
            
            Fqhat=qhat+isnonlinear*[zeros(m,1);compute_fnl(DS,qhat(m+1:end),qdhat(m+1:end))];
            q(:,i)=q(:,i-1)+detT*(qd(:,i-1)+qdhat)/2;
            qd(:,i)=qd(:,i-1)-M\C*(qd(:,i-1)+qdhat)/2*detT-...
                M\K*(Fqxn+Fqhat)/2*detT+M\dF;
                if isnan(q(:,i))
                    error('Forward step does not converge, please decrease the time step')
                end
            %%
%                 dW=zeros(m+n,1); dW(m)=detu;
                dW=dF;
                qhat =(M+1/2*C*detT+1/4*K*detT^2)\(1/2*detT*(M*q2d(:,i-1)+1/2*dW)+(M+1/2*C*detT)*q2(:,i-1));
                qdhat=(M+1/2*C*detT)\(M*q2d(:,i-1)-1/2*K*qhat*detT+1/2*dW);
                q2(:,i)=q2(:,i-1)+detT*qdhat;
                q2d(:,i)=q2d(:,i-1)-M\(C*qdhat+K*qhat)*detT+2*M\dW;
                disp(['time integration completed: ', num2str(i/N*100), '%'])    
%                 disp(['time integration completed: ', num2str(i/N*100), '%'])          
        end
X=q(m+1:end,:); V=qd(m+1:end,:);

X2=q2(m+1:end,:); V2=q2d(m+1:end,:);
[w1,G1]=crossPSDestimator(X(1,:),X(1,:),t);[w2,G2]=crossPSDestimator(X2(1,:),X2(1,:),t);
%%
figure
plot(w1,G1)
hold on
plot(w2,G2)
xlim([0, 20])
legend('heun','implicit')
%%
figure
plot(t,X(1,:))
hold on
plot(t,X2(1,:))
legend('heun','implicit')

%%
function [w,Gxy]=crossPSDestimator(x,y,t)
%%% Reference: Preunont Andre, Random vibration and Spectral Analysis Chp 12.9
%%%% numebr of sampling points 
N = length(x); 
%%%% time period of generation
T0=max(t);
%%%% sampling frequency, time spacing
f0=1/T0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cx=fft(x)/N;Cy=fft(y)/N;Cys=conj(Cy);

Gxy=T0*(Cx.*Cys)/2/pi;
w=(1:length(Cx))*f0*2*pi;

end