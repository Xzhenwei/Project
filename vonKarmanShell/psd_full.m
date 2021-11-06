function psd_full (epsilon)

nDiscretization = 3; % Discretization parameter 30 (#DOFs is proportional to the square of this number)

[M,C,K,fnl,outdof,out] = build_model(nDiscretization);
n = length(M); % number of degrees of freedom

[filterPSD] = build_stochasticF(outdof,n,epsilon);
% Dynamical system setup 
SS = StochasticSystem();   
set(SS,'filterPSD',filterPSD,'linear',false)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS.SSOptions,'ssMethod','indirect','tol',1e-10)
% set(DS.Options,'Emax',5,'Nmax',10,'notation','tensor')
nRealization=2;
T0=30; %% PSD frequency domain resolution is ~ 1/T0
nPoints=2^16; %% control the accuracy of numerical differential equation


clusterRun=false; %% if the script is run on local or cluster.

PSDpair=[out,out];

SS.add_random_forcing(nRealization, T0, nPoints,outdof);

freq_range=[145 155];
[w_linear,linear_analytic] = SS.compute_linear_PSD(PSDpair,freq_range,false);
%%
% *Choose Master subspace (perform resonance analysis)*

S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
set(S.PSDOptions, 'nPointfilter', 8)
masterModes = [1,2];
S.choose_E(masterModes);
order = 5; % Approximation order
%%
S.ssmSEulerTimeDisp = false;
tic
[wss,ssmPSD]=S.extract_PSD(PSDpair, order,'filter',freq_range,clusterRun);
time_ssm=toc;
disp([num2str(time_ssm),' amount of time'])
%%
char=['SSMEp',num2str(epsilon),'.mat'];
save(char);



end