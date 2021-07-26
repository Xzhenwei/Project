%% kappa=1;
n=2; kappa=20;
kappa=0;
[M,C,K,fnl,f_0] = build_model(kappa,n);
% LINEAR QUARTER CAR with Dimension 2;
nRealization=10;
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

clusterRun=false; %% if the script is run on local or cluster.
method="filter ImplicitMidPoint";
PSDpair=[n,n];
%%
[w,outputPSD] = SS.monte_carlo_average(method,PSDpair,nRealization,clusterRun);
%%
% linear_analytic=SS.compute_linear_PSD(SS.input.omega,SS.input.PSD);
%% this section plots the result of the full system response and linear analytical
% figure
% plot(w,outputPSD(1,:),'linewidth',1)
% hold on
% plot(w,linear_analytic(PSDpair(1),:),'linewidth',2)
% hold on
% xlim([0,15])
% xlabel('\Omega frequency')
% ylabel('Power Density')
% legend(strcat('PSD solved by ', ": ",method), 'linear analytic')
% legend('boxoff')
% title('Power Spectral Density of X1')
% grid on
%% SSM setting
S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
masterModes = [1,2];
S.choose_E(masterModes);
order = 5; % SSM approximation order
%%
[wss,ssmPSD]=S.extract_PSD(PSDpair, order,"filter",clusterRun);
%%
figure
plot(wss,ssmPSD,'linewidth',1,'DisplayName','SSM')
hold on
plot(w,outputPSD(1,:),'linewidth',1,'DisplayName','Full System Simulation')
hold on
% plot(w,linear_analytic(PSDpair(1),:),'linewidth',1,'DisplayName','linear analytic')
% xline(firts_res,'-',{'First Resonance'},'linewidth',1.5);
legend
xlim([0,20]);
grid on
