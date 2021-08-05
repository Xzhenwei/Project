function vonKarmanEuler(epsilon)

nElements = 10;
% epsilon = 1e-3;
%% generate model

[M,C,K,fnl,outdof,eMass] = build_model(nElements);
n = length(M);
%%%
nRealization=30;
resol = 3;
T0=resol*100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=resol*2^14; %% control the accuracy of numerical differential equation
[filterPSD, stochastic_f] = build_stochasticF(eMass,n,epsilon);
%% Dynamical system setup  
SS = StochasticSystem();

set(SS,'filterPSD',filterPSD,'linear',false)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl,'gFactor',-eMass);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS.SSOptions,'ssMethod','indirect')
SS.add_random_forcing(nRealization, T0, nPoints,outdof);

%%%%%%%% Above is forcing setting and set to DynamicalSystem class
clusterRun=true; %% if the script is run on local or cluster.
method="filter ImplicitMidPoint";
PSDpair=[n-1,n-1];
%%
[V, D, W] = SS.linear_spectral_analysis();
firts_res=abs(imag(D(1)));
%%
SS.sdeImpTimeDisp = false;
[w,outputPSD] = SS.monte_carlo_average(method,PSDpair,nRealization,clusterRun);
disp('Full system calculation finished')
%% SSM setting
S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
masterModes = [1,2];
S.choose_E(masterModes);
order = 5;
freq_range=[4 7]; % depend on res

S.ssmSEulerTimeDisp = false;
[wss,ssmPSD]=S.extract_PSD(PSDpair, order,'filter heun',freq_range,clusterRun);
[w_linear, linear_analytic]=SS.compute_linear_PSD(PSDpair,freq_range);

char=['vonEpsilon',num2str(epsilon),'.mat'];
save(char,'-mat')

