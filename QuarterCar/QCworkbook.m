%% kappa=1;
n=2; kappa=20;
[M,C,K,fnl,f_0] = build_model(kappa,n);
% LINEAR QUARTER CAR with Dimension 2;
nRealization=1;
T0=100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=2^14; %% control the accuracy of numerical differential equation
epsilon=1000; %% forcing magnitude
[filterPSD,forcingdof,IC,stochastic_f] = build_stochasticF(n,nPoints,epsilon);
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl,'filterPSD',filterPSD,'stochastic_f',stochastic_f,'linear',false); 
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex');
DS.add_random_forcing(nRealization, T0, nPoints,forcingdof);
A=DS.A;B=DS.B;
%%%%%%%% Above is forcing setting and set to DynamicalSystem class

method="filterImplicitMidPoint";
PSDpair=[1,1; 2, 2; 1, 2];
%%
[w,outputPSD] = DS.sde_solver2(method,IC,PSDpair);
%%
linear_analytic=DS.compute_linear_PSD(DS.sde.analytic.omega,DS.sde.analytic.PSD); %% 
%% this section plots the result of the full system response and linear analytical
figure
plot(w,outputPSD(1,:),'linewidth',1)
hold on
plot(w,linear_analytic(1,:),'linewidth',2)
xlim([0,15])
xlabel('\Omega frequency')
ylabel('Power Density')
legend('PSD solved by Implicit Newmark','linear analytic response')
title('Power Spectral Density of X1')
%% this computes the spectral properties
[V, D, W] = DS.linear_spectral_analysis();
firts_res=abs(imag(D(1)));
%% SSM setting
S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
masterModes = [1,2];
S.choose_E(masterModes);
order = 5; % SSM approximation order
[W0, R0] = S.compute_whisker(order);
ssmPSD=S.extract_PSD([0 20], order,"filter");
%%
figure
plot(ssmPSD.omega,ssmPSD.PSD,'linewidth',1)
hold on
plot(w,outputPSD(1,:),'linewidth',1)
% hold on
% plot(w,linear_analytic(1,:),'linewidth',1)
% legend('SSM','nonlinear simulation','linear analytic')
xlim([0,20]);
grid on
