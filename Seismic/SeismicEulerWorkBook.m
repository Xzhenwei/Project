n = 10;
m = 3;
k = 2000;
c = 50;
kappa2 = 0;
kappa3 = 10;

[M,C,K,fnl] = build_model(n,m,c,k,kappa2,kappa3);
nRealization=20;
T0=100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=2^14; %% control the accuracy of numerical differential equation
epsilon=1000; %% forcing magnitude
[filterPSD,forcingdof,stochastic_f] = build_stochasticF(nPoints,epsilon);
%%%%%%%%
DS = DynamicalSystem();

SS = StochasticSystem();

set(SS,'filterPSD',filterPSD,'linear',false)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl,'n',n);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
SS.add_random_forcing(nRealization, T0, nPoints,forcingdof);

%%%%%%%% Above is forcing setting and set to DynamicalSystem class

method="filter ImplicitMidPoint";
PSDpair=[n,n];
%%
[V, D, W] = SS.linear_spectral_analysis();
clusterRun=true;
%%
tic
[w,outputPSD] = SS.monte_carlo_average(method,PSDpair,nRealization,clusterRun);
time_sde=toc;
disp(['Total number of ',num2str(nRealization),'# takes ',...
    num2str(time_sde),' amount of time'])
%%
linear_analytic=SS.compute_linear_PSD(SS.input.omega,SS.input.PSD);

%% SSM setting
S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
masterModes = [1,2];
S.choose_E(masterModes);
order = 5; % SSM approximation order

[wss,ssmPSD]=S.compute_ssmPSD(PSDpair, order,"filter",clusterRun);
