function linear_and_full (epsilon)

tol = 1e-12 ;
nDiscretization = 10; % Discretization parameter 30 (#DOFs is proportional to the square of this number)
nRealization=24;
T0=30; %% PSD frequency domain resolution is ~ 1/T0
nPoints=2^16; %% control the accuracy of numerical differential equation




[M,C,K,fnl,outdof,out] = build_model(nDiscretization);
n = length(M); % number of degrees of freedom

fnl_l = {sptensor([n,n,n]),sptensor([n,n,n,n])};
[filterPSD] = build_stochasticF(outdof,n,epsilon);
% Dynamical system setup
SS = StochasticSystem();
set(SS,'filterPSD',filterPSD,'linear',false)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS.SSOptions,'ssMethod','indirect','tol',tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SS_l = StochasticSystem();
set(SS_l,'filterPSD',filterPSD,'linear',true)
set(SS_l,'M',M,'C',C,'K',K,'fnl',fnl_l);
set(SS_l.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS_l.SSOptions,'ssMethod','indirect','tol',tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clusterRun=true; %% if the script is run on local or cluster.
method="filter ImplicitMidPoint";
PSDpair=[out,out];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SS.add_random_forcing(nRealization, T0, nPoints,outdof);
SS_l.add_random_forcing(nRealization, T0, nPoints,outdof);
[V,D,W] = SS.linear_spectral_analysis();
freq_range=[140 170];
SS.sdeImpTimeDisp = false; SS_l.sdeImpTimeDisp = false;

tic
[w_l,outputPSD_l] = SS_l.monte_carlo_average(method,PSDpair,nRealization,clusterRun);
time_sde=toc;
disp(['Total number of ',num2str(nRealization),'# realization takes ',...
    num2str(time_sde),' amount of time in linear simulation'])


tic
[w,outputPSD] = SS.monte_carlo_average(method,PSDpair,nRealization,clusterRun);
time_sde=toc;
disp(['Total number of ',num2str(nRealization),'# realization takes ',...
    num2str(time_sde),' amount of time in full nonlinear simulation'])
S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
set(S.PSDOptions, 'nPointfilter', 8)
masterModes = [1,2];
S.choose_E(masterModes);
order = 5; % Approximation order
%%
char=['LinAndFullEp',num2str(epsilon),'T',num2str(T0),'nP',num2str(log(nPoints))];
save(char);
%%
S.ssmSEulerTimeDisp = false;
tic
[wss,ssmPSD]=S.extract_PSD(PSDpair, order,'filter',freq_range,clusterRun);
time_ssm=toc;
disp([num2str(time_ssm),' amount of time'])

char=['LinAndFullEp',num2str(epsilon),'T',num2str(T0),'nP',num2str(log(nPoints)/log(2)),'.mat'];
save(char);


end
