function linear_and_full (epsilon)

clear all; close all; clc

nDiscretization = 10; % Discretization parameter 30 (#DOFs is proportional to the square of this number)
epsilon = 0.1;

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

nRealization=1;
T0=100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=2^14; %% control the accuracy of numerical differential equation


clusterRun=true; %% if the script is run on local or cluster.
method="filter ImplicitMidPoint";
PSDpair=[out,out];

SS.add_random_forcing(nRealization, T0, nPoints,outdof);


freq_range=[145 155];
tic
[w_linear,linear_analytic] = SS.compute_linear_PSD(PSDpair,freq_range,clusterRun);
timeconsume = toc;
char=['linearAnalyticEp',num2str(epsilon),'.mat'];
save(char,'w_linear','linear_analytic','timeconsume');



end