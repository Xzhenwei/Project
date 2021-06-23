function OUTPUT = extract_PSD(obj, parRange, ORDER, method)
%  EXTRACT_FRC This function extracts the power spectral density (PSD) for
%  systems under stochastic forcing. The PSD computation is based
%  on SSM computation. An appropriate SSM is constructed based on the
%  resonant spectrum of system. The response is computed for the reduced system firstly. 
%  The obtained response is finally mapped back to physical coordinates and 
%  compute the PSD. 
%  
%
% parRange: range of frequency
% order:    order of SSM expansion to be used for FRC computation

% order=[3 5];
% for j = 1:numel(ORDER)
%     order = ORDER(j);  % SSM approximation order
order = ORDER;
[W0, R0] = obj.SSM.compute_whisker(order);
%% computation of the reduced dynamics
m = obj.SSM.dimManifold; 

% the resolution of solving 2-dim system is higher than original system in
% order to capture the high frequency information
num_points=obj.nPoints*2^2;
T=obj.timeSpan;
p0 = zeros(m,1);

detT=T/num_points; t=0:detT:T; 

Wnode=obj.SSM.E.adjointBasis';

MontCarlo=obj.nRealization;   Gzz=0;
n=obj.System.n;


%% backward euler
for l=1:MontCarlo
% S.E.basis/adjointBasis
%     for i=1:num_points % implicit Euler
%         maxiter=100;
%         options = optimoptions('fsolve','Display','off');
%         pdot=expand_coefficients(DS,R0,Wnode,m,t(i), p(:,i));
%         pt=p(:,i)+pdot*detT;
%         f_p=@(x) expand_coefficients(DS,R0,Wnode,m,t(i+1), x)*detT+p(:,i)-x;
%         p(:,i+1)=fsolve(f_p,pt, options);
% 
%         disp(['time integration completed: ', num2str(i/num_points*100), '%']);
%     end
% z=zeros(2*n,length(p));% num_points+1
%     for i=1:length(p)
%         z(:,i) = expand_autonomous(W0,2*n, p(:,i));
%     end
% [w,Gz]=crossPSDestimator(z(1,:),z(1,:),t);
% Gzz=Gzz+Gz;
%% indirect filter method
switch lower(method)
    case "filter"
        PSD=obj.filterPSD;
        Mz=PSD.Mz;
        PSD.G=zeros(n,length(Mz)); %G(end,:)=ones(1,length(Mz));
        f=length(Mz);
        p=obj.indirect_Euler_SSM(num_points,T,PSD,f,m,Wnode,R0);

        z=zeros(2*n,length(p));
        for i=1:length(p)
            z(:,i) = expand_autonomous(W0,2*n, p(:,i));
        end
        [w,Gz]=crossPSDestimator(z(1,:),z(1,:),t);
        Gzz=Gzz+Gz;
    
    otherwise
%%
        obj.System.Fsto=obj.System.generate_stochastic('indirect');

%     %% ODE45 solver
        f_p=@(t,p) expand_coefficients(obj.System,R0,Wnode,m,t, p);
        opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
        [t_45,p_45] = ode45(f_p,t,p0, opts);
        p_45=p_45';
        z=zeros(2*n,length(p_45));
        for i=1:length(p_45)
            z(:,i) = expand_autonomous(W0,2*n, p_45(:,i));
        end
        [w,Gz]=crossPSDestimator(z(1,:),z(1,:),t_45);
        Gzz=Gzz+Gz;
end
%% calculate full system PSD       

disp(['reamaining number of Monte Carlo simulation: ', num2str((MontCarlo-l))])
end

Gzz=Gzz/MontCarlo;

%% calculating non-autonomous analytically
        omega=obj.input.omega;
        forcePSD=obj.input.PSD;
    
    B=obj.System.B; M=obj.System.M; C=obj.System.C; K=obj.System.K; 
    F_psd=zeros(n,n);
    G=eye(2*n)-B*obj.E.basis*Wnode;
    G11=G(1:n,1:n);
    Z11=zeros(1,length(omega)); 
    forcingdof=obj.forcingdof;
for j=1:length(omega)
% second order system
    F_psd(forcingdof,forcingdof)=forcePSD(j);
    Hw=inv(-omega(j)^2*M+1i*omega(j)*C+K);
    Z_full=(-omega(j)^2*M+1i*omega(j)*C+K)\G11*F_psd*G11'*Hw';
    Z11(j)=norm(Z_full(1,1));

end
%%
num_points=obj.num_points;
%%%Adding PSD
We=zeros(num_points+1,1)'; 

for i=3:num_points+1
    We(i)= interp1(w,Gzz,omega(i)) +Z11(i);
end
OUTPUT.PSD=We;
OUTPUT.omega=omega;
%%
figure
plot(omega,We,'linewidth', 1.5)
xlim(parRange);
grid on

%%

function S = expand_coefficients(obj,S0,U,n, t, p)
% m = numel(p);
S = zeros(n,1);

% expand autonomous coefficients
for k = 1:length(S0)
    S =  S + expand_multiindex(S0{k},p);
end

% add nonautonomous contribution
    S = S + U * [compute_fstochastic(obj,t); zeros(obj.n,1)]; 
end

function S = expand_autonomous(W0,n, p)
S = zeros(n,1);
% expand autonomous coefficients
for k = 1:length(W0)
    S =  S + real(expand_multiindex(W0{k},p));
end

end

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



end