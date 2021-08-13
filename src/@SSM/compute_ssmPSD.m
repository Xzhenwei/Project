function [w,Gzz] = compute_ssmPSD(obj, PSDpair, W0, R0, method)
%  EXTRACT_FRC This function extracts the power spectral density (PSD) for
%  systems under stochastic forcing. The PSD computation is based
%  on SSM computation. An appropriate SSM is constructed based on the
%  resonant spectrum of system. The response is computed for the reduced system firstly. 
%  The obtained response is finally mapped back to physical coordinates and 
%  compute the PSD.


%% computation of the reduced dynamics
m = obj.dimManifold; 
%%% 
% the resolution of solving 2-dim system is higher than original system in
% order to capture the high frequency information
num_points = obj.System.nPoints*2^2; % 3 when filter was used.

T=obj.System.timeSpan;
p0 = sparse(m,1);

detT=T/num_points; t=0:detT:T; 

Wnode=obj.E.adjointBasis';

n=obj.System.n;

nOutput=size(PSDpair,1);
%%
Gzz = sparse(nOutput,num_points+1);
for j=1:nOutput

%% indirect filter method
switch lower(method)
    case 'filter'
        PSD=obj.System.filterPSD;
        Mz=PSD.Mz;
        f=length(Mz);
        p=ssm_Implicit_solver(obj, num_points,T,PSD,f,m,Wnode,R0);
%         p=ssm_Heun_solver(obj, num_points,T,PSD,f,m,Wnode,R0);
        
%         z=sparse(2*n,length(p));
        z1=zeros(1,length(p)); z2=z1;
        for i=1:length(p)
            AS = expand_autonomous(W0,2*n, p(:,i));
            z1(i) = AS(PSDpair(j,1));
            z2(i) = AS(PSDpair(j,2));
        end
        [w,Gz]=crossPSDestimator(z1,z2,t);
        Gzz(j,:)=Gz;
     case 'filter heun'
        PSD=obj.System.filterPSD;
        Mz=PSD.Mz;
        f=length(Mz);
        p=ssm_Heun_solver(obj, num_points,T,PSD,f,m,Wnode,R0);
        
        z1=zeros(1,length(p)); z2=z1;
        for i=1:length(p)
            AS = expand_autonomous(W0,2*n, p(:,i));
            z1(i) = AS(PSDpair(j,1));
            z2(i) = AS(PSDpair(j,2));
        end
        [w,Gz]=crossPSDestimator(z1,z2,t);
        Gzz(j,:)=Gz;
    
    otherwise
%%
        obj.System.Fsto=obj.System.generate_stochastic();

%     %% ODE45 solver
%         detT=T/num_points*4; t=0:detT:T;
        f_p=@(t,p) expand_DE(obj,R0,Wnode,m,t, p,T);
        opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
        t_span = t;
        [t_45,p_45] = ode45(f_p,t_span ,p0, opts);
        p_45=p_45';
        
        z1=zeros(1,length(p_45)); z2=z1;
        for i=1:length(p_45)
            AS = expand_autonomous(W0,2*n, p_45(:,i));
            z1(i) = AS(PSDpair(j,1));
            z2(i) = AS(PSDpair(j,2));
        end
        [w,Gz]=crossPSDestimator(z1,z2,t_45);
        Gzz(j,:)=Gz;
        
end
end



%%

function S = expand_DE(obj,S0,U,n, t, p,T)
% m = numel(p);
S = sparse(n,1);

% expand autonomous coefficients
    for k = 1:length(S0)
        S =  S + expand_multiindex(S0{k},p);
    end
% add nonautonomous contribution
    S = S + U * [obj.System.compute_fstochastic(t); sparse(obj.System.n,1)];
    
    if obj.ssmSEulerTimeDisp
    disp(['Integration completed: ', num2str(t/T*100),...
        '%'])
    end
    
end

function S = expand_autonomous(W0,n, p)
S = sparse(n,1);
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