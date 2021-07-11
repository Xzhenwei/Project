
% *system parameters*

nElements = 10;
epsilon = 0;
%% generate model

[M,C,K,fnl,f_0,outdof] = build_model(nElements);
n = length(M);
%%%
nRealization=10;
T0=100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=2^14; %% control the accuracy of numerical differential equation
% epsilon=1000; %% forcing magnitude
[filterPSD,forcingdof,stochastic_f] = build_stochasticF(nPoints,epsilon);
%% Dynamical system setup 

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
tic
[w,outputPSD] = SS.monte_carlo_average(method,PSDpair,nRealization);
time_sde=toc;
disp(['Total number of ',num2str(nRealization),'# takes ',...
    num2str(time_sde),' amount of time'])

%% Linear Modal analysis and SSM setup

[V,D,W] = SS.linear_spectral_analysis();
linear_analytic=SS.compute_linear_PSD(SS.input.omega,SS.input.PSD);

