%% kappa=1;
n=2; kappa=200;
[M,C,K,fnl,f_0] = build_model(kappa,n);
% LINEAR QUARTER CAR with Dimension 2;
nRealization=1;
T0=100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=2^14; %% control the accuracy of numerical differential equation
epsilon=1000; %% forcing magnitude
[filterPSD,forcingdof,IC,stochastic_f] = build_stochasticF(n,nPoints,epsilon);
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex');

%%
SDE=StochasticSystem(DS,'System');
set(SDE,'filterPSD',filterPSD,'linear',false)
SDE.add_random_forcing(nRealization, T0, nPoints,forcingdof);

%%%%%%%% Above is forcing setting and set to DynamicalSystem class

method="filter ImplicitMidPoint";
PSDpair=[1,1; 2, 2; 1, 2];
%%
[w,outputPSD] = SDE.sde_solver(method,PSDpair);
%%
linear_analytic=SDE.compute_linear_PSD(SDE.input.omega,SDE.input.PSD); 
%% this section plots the result of the full system response and linear analytical
figure
plot(w,outputPSD(1,:),'linewidth',1)
hold on
plot(w,linear_analytic(1,:),'linewidth',2)
hold on
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
ssmPSD=SDE.extract_PSD([0 20], order,"");
%%
figure
plot(ssmPSD.omega,ssmPSD.PSD,'linewidth',1)
hold on
plot(w,outputPSD(1,:),'linewidth',1)
hold on
plot(w,linear_analytic(1,:),'linewidth',1)
legend('SSM','nonlinear simulation','linear analytic')
xlim([0,20]);
grid on
