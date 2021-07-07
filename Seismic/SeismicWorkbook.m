n = 10;
m = 3;
k = 2000;
c = 50;
kappa2 = 0;
kappa3 = 0;

[M,C,K,fnl] = build_model(n,m,c,k,kappa2,kappa3);
nRealization=20;
T0=100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=2^14; %% control the accuracy of numerical differential equation
epsilon=1000; %% forcing magnitude
[filterPSD,forcingdof,stochastic_f] = build_stochasticF(nPoints,epsilon);
%%%%%%%%
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')

SDE=StochasticSystem(DS,'System');
set(SDE,'filterPSD',filterPSD,'linear',false)
SDE.add_random_forcing(nRealization, T0, nPoints,forcingdof);
%%%%%%%% Above is forcing setting and set to DynamicalSystem class

method="filter ImplicitMidPoint";
PSDpair=[n,n];
%%
[V, D, W] = DS.linear_spectral_analysis();
firts_res=abs(imag(D(1)));
%%
[w,outputPSD] = SDE.sde_solver(method,PSDpair);
%%
linear_analytic=SDE.compute_linear_PSD(SDE.input.omega,SDE.input.PSD);
%% this section plots the result of the full system response and linear analytical
figure
plot(w,outputPSD(1,:),'linewidth',1)
hold on
plot(w,linear_analytic(n,:),'linewidth',2)
xlim([0,15])
xlabel('\Omega frequency')
ylabel('Power Density')
legend(strcat('PSD solved by ', ": ",method), 'linear analytic')
legend('boxoff')
title('Power Spectral Density of X1')
grid on

%% SSM setting
S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
masterModes = [1,2];
S.choose_E(masterModes);
order = 5; % SSM approximation order
%%
SDE.SSM=S;
ssmPSD=SDE.extract_PSD([0 20], order,"filter");
%%
figure
plot(ssmPSD.omega,ssmPSD.PSD,'linewidth',1,'DisplayName','SSM')
hold on
plot(w,outputPSD(1,:),'linewidth',1,'DisplayName','Full System Simulation')
% hold on
% plot(w,linear_analytic(n,:),'linewidth',1,'DisplayName','linear analytic')
% xline(firts_res,'-',{'First Resonance'},'linewidth',1.5);
legend
xlim([0,10]);
grid on
