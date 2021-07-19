
nElements = 10;
epsilon = 1e-2;

[M,C,K,fnl,outdof] = build_model(nElements);
n = length(M);

nRealization=20;
T0=100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=2^14; %% control the accuracy of numerical differential equation
% epsilon=1000; %% forcing magnitude
[filterPSD, stochastic_f] = build_stochasticF(nPoints,epsilon);
%% Dynamical system setup  
SS = StochasticSystem();

set(SS,'filterPSD',filterPSD,'linear',true)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS.SSOptions,'ssMethod','indirect')
SS.add_random_forcing(nRealization, T0, nPoints,outdof);

%%%%%%%% Above is forcing setting and set to DynamicalSystem class
clusterRun=true; %% if the script is run on local or cluster.
method="filter ImplicitMidPoint";
PSDpair=[n,n];
%%
[V, D, W] = SS.linear_spectral_analysis();
%%
SS.sdeImpTimeDisp = false;
tic
[w,outputPSD] = SS.monte_carlo_average(method,PSDpair,nRealization,clusterRun);
time_sde=toc;
disp(['Total number of ',num2str(nRealization),'# realization takes ',...
    num2str(time_sde),' amount of time'])

%% SSM setting
S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
masterModes = [1,2];
S.choose_E(masterModes);
order = 5; % SSM approximation order
%%
S.ssmSEulerTimeDisp=false;
[wss,ssmPSD]=S.compute_ssmPSD(PSDpair, order,"",clusterRun);
%%
linear_analytic=SS.compute_linear_PSD(SS.input.omega,SS.input.PSD);
