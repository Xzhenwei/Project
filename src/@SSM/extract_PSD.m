function OUTPUT = extract_PSD(obj, PSDpair, ORDER, method,clusterRun)
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

%% computation of the reduced dynamics
n = obj.System.n;

[w,Gzz] = obj.compute_ssmPSD( PSDpair, ORDER, method);

Wnode = obj.E.adjointBasis';

assert(obj.System.nRealization > 1, ' Please use compute_ssmPSD for one simulation')

if clusterRun
    euler = parcluster('local');
    pool = parpool(euler,24);
else
    pool = parpool('local',2);
end

parfor i=1:obj.System.nRealization-1
    [~,outputPSD] = obj.compute_ssmPSD( PSDpair, ORDER, method);
    Gzz = Gzz+outputPSD;
end
Gzz = Gzz/obj.System.nRealization;

pool.delete()


%% calculating non-autonomous analytically
        omega=obj.System.input.omega;
        forcePSD=obj.System.input.PSD;
    
    B=obj.System.B; M=obj.System.M; C=obj.System.C; K=obj.System.K; 
    F_psd=zeros(n,n);
    G=eye(2*n)-B*obj.E.basis*Wnode;
    G11=G(1:n,1:n);
    Z=zeros(1,length(omega)); 
    forcingdof=obj.System.forcingdof;
for j=1:length(omega)
% second order system
    F_psd(forcingdof,forcingdof)=forcePSD(j);
    Hw=inv(-omega(j)^2*M+1i*omega(j)*C+K);
    Z_full=(-omega(j)^2*M+1i*omega(j)*C+K)\G11*F_psd*G11'*Hw';
    Z(j)=norm(Z_full(PSDpair(1),PSDpair(2)));

end
%%
nPoints=obj.System.nPoints;
%%%Adding PSD
We=zeros(nPoints+1,1)'; 

for i=3:nPoints+1
    We(i)= interp1(w,Gzz,omega(i)) +Z(i);
end
OUTPUT.PSD=We;
OUTPUT.omega=omega;


end