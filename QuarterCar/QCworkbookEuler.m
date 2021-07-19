%% kappa=1;
n=2; kappa=200;
kappa=0;
[M,C,K,fnl,f_0] = build_model(kappa,n);
% LINEAR QUARTER CAR with Dimension 2;
nRealization=20;
T0=100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=2^13; %% control the accuracy of numerical differential equation
epsilon=1000; %% forcing magnitude
[filterPSD,forcingdof,IC,stochastic_f] = build_stochasticF(n,nPoints,epsilon);
DS = DynamicalSystem();

SS = StochasticSystem();

set(SS,'filterPSD',filterPSD,'linear',false)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
SS.add_random_forcing(nRealization, T0, nPoints,forcingdof);

%%%%%%%% Above is forcing setting and set to DynamicalSystem class

clusterRun=true; %% if the script is run on local or cluster.
method="filter ImplicitMidPoint";
PSDpair=[n,n];
%%
[w,outputPSD] = SS.sde_solver(method,PSDpair);
%%
linear_analytic=SS.compute_linear_PSD(SS.input.omega,SS.input.PSD);
%% SSM setting
S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
masterModes = [1,2];
S.choose_E(masterModes);
order = 5; % SSM approximation order
%%
[wss,ssmPSD]=S.compute_ssmPSD(PSDpair, order,"filter",clusterRun);
save('QCworkspace.mat')
