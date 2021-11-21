function vonKarmanEuler(epsilon)

% *system parameters*
nElements = 10;
%% generate model
[M,C,K,fnl,outdof,eMass] = build_model(nElements);
n = length(M);
%%%
nRealization = 48 ;
T0= 5*100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=5*2^14; %% control the accuracy of numerical differential equation

[filterPSD, stochastic_f] = build_stochasticF(eMass,n,epsilon);
%% Dynamical system setup  
SS = StochasticSystem();

set(SS,'filterPSD',filterPSD,'linear',false)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl,'gFactor',-eMass);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS.SSOptions,'ssMethod','indirect')
SS.add_random_forcing(nRealization, T0, nPoints,outdof);
SS.inputForcingType = 'disp';
%%%%%%%% Above is forcing setting and set to DynamicalSystem class
clusterRun=true; %% if the script is run on local or cluster.
method="filter ImplicitMidPoint";
PSDpair=[n-1,n-1];
%%
[V, D, W] = SS.linear_spectral_analysis();
first_res=abs(imag(D(1)));
%%
tic
[w,outputPSD] = SS.monte_carlo_average(method,PSDpair,nRealization,clusterRun);
time_sde=toc;
disp(['Total number of ',num2str(nRealization),'# realization takes ',...
    num2str(time_sde),' amount of time'])
%%
tic
[w_galerkin, PSD_galerkin] =  SS.monte_carlo_galerkin(method, PSDpair);
disp(['Total number of ',num2str(nRealization),'# realization in Galerkin projection takes ',...
    num2str(time_sde),' amount of time'])
%% SSM setting
S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
set(S.PSDOptions, 'nPointfilter', 1)
masterModes = [1,2];
S.choose_E(masterModes);
order = 5; % SSM approximation order
%%
freq_range=[3 10];
S.ssmSEulerTimeDisp = true;
tic
[wss,ssmPSD]=S.extract_PSD(PSDpair, order,'filter heun',freq_range,clusterRun);
time_ssm=toc;
disp([num2str(time_ssm),'s amount of time'])

%%
[w_linear, linear_analytic] = SS.compute_linear_PSD(PSDpair,freq_range);

%%
char=['vonBeamEpsilon',num2str(epsilon),'.mat'];
save(char,'-mat')

end