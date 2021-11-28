function ShellEuler(epsilon)

clear all; close all; clc

nDiscretization = 10; % Discretization parameter 30 (#DOFs is proportional to the square of this number)
% epsilon = 0.02;
%% generate model

[M,C,K,fnl,outdof,out] = build_model(nDiscretization);
n = length(M); % number of degrees of freedom
disp(['Number of degrees of freedom = ' num2str(n)])
disp(['Phase space dimensionality = ' num2str(2*n)])
%%
[filterPSD] = build_stochasticF(outdof,n,epsilon);

SS = StochasticSystem();   
set(SS,'filterPSD',filterPSD,'linear',false)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS.SSOptions,'ssMethod','indirect','tol',1e-12)
% set(DS.Options,'Emax',5,'Nmax',10,'notation','tensor')
nRealization = 24;
T0 = 30; %% PSD frequency domain resolution is ~ 1/T0
nPoints = 2^16; %% control the accuracy of numerical differential equation
%% 
% We assume periodic forcing of the form
% 
% $$\mathbf{f}^{ext}(\phi) = \mathbf{f}_0\cos(\phi)=\frac{\mathbf{f}_0}{2}e^{i\phi} 
% + \frac{\mathbf{f}_0}{2}e^{-i\phi}  $$
% 
%%%%%%%% Above is forcing setting and set to DynamicalSystem class
clusterRun = true; %% if the script is run on local or cluster.
method = "filter ImplicitMidPoint";
PSDpair = [out,out];

SS.add_random_forcing(nRealization, T0, nPoints,outdof);

%% Linear Modal analysis and SDE
[V,D,W] = SS.linear_spectral_analysis();
firts_res = abs(imag(D(1)));
SS.sdeImpTimeDisp = false;
tic
[w,outputPSD] = SS.monte_carlo_average(method,PSDpair,nRealization,clusterRun);
time_sde = toc;
disp(['Total number of ',num2str(nRealization),'# realization takes ',...
    num2str(time_sde),' amount of time'])
%% Linear Analytic
freq_range = [30 55];
[w_linear,linear_analytic] = SS.compute_linear_PSD(PSDpair,freq_range);

%%
% *Choose Master subspace (perform resonance analysis)*
SS.nPoints = 2^17;
S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
set(S.PSDOptions, 'nPointfilter', 16)
masterModes = [1,2];
S.choose_E(masterModes);
order = 5; % Approximation order
%%
S.ssmSEulerTimeDisp = false;
tic
[wss,ssmPSD] = S.extract_PSD(PSDpair, order,'filter heun',freq_range,clusterRun);
time_ssm = toc;
disp([num2str(time_ssm),' amount of time'])

%%
[U, OmegaSq, NOT_CONVERGED] = eigs(sparse(K),sparse(M),2,'smallestabs');
U = U(:,1) ;
% record time and place in table
euler = parcluster('local');
pool = parpool(euler);
PSD_galerkin = 0; w_galerkin = 0;
tic
parfor i=1:nRealization
    
    [w_g,Gz_g] = galerkin_proj(SS, U, filterPSD, method, PSDpair);
    PSD_galerkin = PSD_galerkin+Gz_g;
    w_galerkin = w_galerkin + w_g;
    disp(['number of realizations left:', num2str(nRealization-i)])
    
end
time_glk = toc;
PSD_galerkin = PSD_galerkin/nRealization;
w_galerkin = w_galerkin/nRealization;
pool.delete()

char=[num2str(nDiscretization),'ShellEpsilon',num2str(epsilon),'T',num2str(T0),'nP',num2str(log(nPoints)/log(2)),'.mat'];
save(char)
% ENHANCED MONTE CARLO/ Table of computational times for each examples.
% write a master script

% epsilon 0.1-0.25s

% figure
% plot(w1,outputPSD1,'Displayname','imp')
% hold on
% plot(w2,outputPSD2,'Displayname','heun')
% xlim([120 160])
% legend