function psd_full (epsilon)


nDiscretization = 10; % Discretization parameter 30 (#DOFs is proportional to the square of this number)

[M,C,K,fnl,outdof,out] = build_model(nDiscretization);
n = length(M); % number of degrees of freedom
disp(['Number of degrees of freedom = ' num2str(n)])
disp(['Phase space dimensionality = ' num2str(2*n)])

[filterPSD, stochastic_f] = build_stochasticF(outdof,n,epsilon);
% Dynamical system setup 
SS = StochasticSystem();   
set(SS,'filterPSD',filterPSD,'linear',false)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS.SSOptions,'ssMethod','indirect')

nRealization=30;
T0=10; %% PSD frequency domain resolution is ~ 1/T0
nPoints=2^11; %% control the accuracy of numerical differential equation


clusterRun=true; %% if the script is run on local or cluster.
method="filter ImplicitMidPoint";
PSDpair=[out,out];

SS.add_random_forcing(nRealization, T0, nPoints,outdof);


tic
[w,outputPSD] = SS.monte_carlo_average(method,PSDpair,nRealization,clusterRun);
time_sde=toc;
disp(['Total number of ',num2str(nRealization),'# realization takes ',...
    num2str(time_sde),' amount of time'])
char=['FullEp',num2str(epsilon),'.mat'];
save(char,'w','outputPSD','time_sde');



end