n = 10;
m = 7;
k = 4555;
c = 90;
kappa2 = 0;
kappa3 = 1.5;

[M,C,K,fnl] = build_model(n,m,c,k,kappa2,kappa3);
nRealization=1;
T0=100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=2^14; %% control the accuracy of numerical differential equation
epsilon=0.01; %% forcing magnitude
[filterPSD,forcingdof,stochastic_f] = build_stochasticF(m,epsilon,n);
%%%%%%%%
SS = StochasticSystem();

set(SS,'filterPSD',filterPSD,'linear',false,'inputForcingType','disp')
set(SS,'M',M,'C',C,'K',K,'fnl',fnl,'gFactor',-m);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS.SSOptions,'ssMethod','indirect')
SS.add_random_forcing(nRealization, T0, nPoints,forcingdof);

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
%%
[w_linear, linear_analytic]=SS.compute_linear_PSD(PSDpair,freq_range);
%% SSM setting
S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
set(S.PSDOptions, 'nPointfilter', 1)
masterModes = [1,2];
S.choose_E(masterModes);
order = 5; % SSM approximation order
%%
tic
[wss,ssmPSD] = S.extract_PSD(PSDpair, order,'filter heun',freq_range,clusterRun);
time_ssm = toc;
disp([num2str(time_ssm),' amount of time'])

%%
clear omega
clear Gxx
omega.w = w; Gxx.Gfull = outputPSD;%outputPSD;
omega.linear = w_linear; Gxx.linear_analytic = linear_analytic;
omega.wss = wss; Gxx.Gss = ssmPSD;
omega.w_galerkin = w_galerkin; Gxx.galerkin = PSD_galerkin;

plot_log_PSD(omega,Gxx,order,PSDpair,[0 7],false)
