function [w,Gzz] = compute_ssmPSD(obj, PSDpair, ORDER, method)

m = obj.dimManifold; 
order = ORDER;
[W0, R0] = obj.compute_whisker(order);

% the resolution of solving 2-dim system is higher than original system in
% order to capture the high frequency information
num_points = obj.System.nPoints*2^2;
T = obj.System.timeSpan;
p0 = zeros(m,1);

detT = T/num_points; t=0:detT:T; 

Wnode = obj.E.adjointBasis';

n = obj.System.n;

switch lower(method)
    case "filter"
        PSD=obj.System.filterPSD;
        Mz=PSD.Mz;
        PSD.G=zeros(n,length(Mz)); %G(end,:)=ones(1,length(Mz));
        f=length(Mz);
        p=obj.indirect_Euler_SSM(num_points,T,PSD,f,m,Wnode,R0);

        z=zeros(2*n,length(p));
        for i=1:length(p)
            z(:,i) = expand_autonomous(W0,2*n, p(:,i));
        end
        [w,Gzz]=crossPSDestimator(z(PSDpair(1),:),z(PSDpair(2),:),t);
    
    otherwise

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
        [w,Gzz]=crossPSDestimator(z(PSDpair(1),:),z(PSDpair(2),:),t_45);
end


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