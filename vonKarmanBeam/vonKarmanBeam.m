
% *system parameters*

nElements = 10;
epsilon = 1;
%% generate model

[M,C,K,fnl,outdof] = build_model(nElements);
n = length(M);
%%%
nRealization=1;
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
clusterRun=false; %% if the script is run on local or cluster.
method="filter ImplicitMidPoint";
PSDpair=[n,n];
%%
[V, D, W] = SS.linear_spectral_analysis();
firts_res=abs(imag(D(1)));
%%
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
[wss,ssmPSD]=S.compute_ssmPSD(PSDpair, order,"filter",clusterRun);
%%
linear_analytic=SS.compute_linear_PSD(SS.input.omega,SS.input.PSD);
%%
figure
plot(wss,ssmPSD,'linewidth',1,'DisplayName','SSM')
hold on
plot(w,outputPSD(1,:),'linewidth',1,'DisplayName','Full System Simulation')
hold on
plot(w,linear_analytic(n,:),'linewidth',1,'DisplayName','linear analytic')
% % xline(firts_res,'-',{'First Resonance'},'linewidth',1.5);
legend
xlim([0,10]);
grid on
