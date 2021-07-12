n = 10;
m = 3;
k = 2000;
c = 50;
kappa2 = 0;
kappa3 = 0;

[M,C,K,fnl] = build_model(n,m,c,k,kappa2,kappa3);
nRealization=2;
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
firts_res=abs(imag(D(1)));
%%
tic
[w,outputPSD] = SS.monte_carlo_average(method,PSDpair,nRealization);
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

ssmPSD=S.compute_ssmPSD(PSDpair, order,"filter");

%% this section plots the result
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
