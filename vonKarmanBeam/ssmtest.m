function ssmtest(epsilon)
nElements = 10;
[M,C,K,fnl,outdof,eMass] = build_model(nElements);
n = length(M);
nRealization = 24 ;
T0= 100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=2^15; %% control the accuracy of numerical differential equation

[filterPSD, ~] = build_stochasticF(eMass,n,epsilon);
%% Dynamical system setup  
SS = StochasticSystem();

set(SS,'filterPSD',filterPSD,'linear',false)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl,'gFactor',-eMass);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS.SSOptions,'ssMethod','indirect')
SS.add_random_forcing(nRealization, T0, nPoints,outdof);
SS.inputForcingType = 'disp';
%%%%%%%% Above is forcing setting and set to DynamicalSystem class
PSDpair=[n-1,n-1];
%% SSM setting
S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
set(S.PSDOptions, 'nPointfilter', 1)
masterModes = [1,2];
S.choose_E(masterModes);
order = 5; % SSM approximation order
%%
freq_range=[3 10];
S.ssmSEulerTimeDisp = false;
tic
[wss,ssmPSD]=S.extract_PSD(PSDpair, order,'filter',freq_range);
time_ssm=toc;
disp([num2str(time_ssm),'s amount of time'])

%%
[w_linear, linear_analytic] = SS.compute_linear_PSD(PSDpair,freq_range);
%%
char=['ssmBeamEpsilon',num2str(epsilon),'.mat'];
save(char,'-mat')
end