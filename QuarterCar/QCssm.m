function QCssm(epsilon)
%% kappa=1;
n=2; kappa=2e4;
% kappa=0;
[M,C,K,fnl,f_0] = build_model(kappa,n);
% LINEAR QUARTER CAR with Dimension 2;
nRealization=48;
T0=50; %% PSD frequency domain resolution is ~ w0= 2*pi* 1/T0
nPoints=2^14; %% control the accuracy of numerical differential equation
[samplePSD,forcingdof,IC,stochastic_f] = build_stochasticF(n,nPoints,epsilon);
DS = DynamicalSystem();

SS = StochasticSystem();

set(SS,'samplePSD',samplePSD,'linear',false)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS.SSOptions,'ssMethod','direct','tol',1e-12)
SS.add_random_forcing(nRealization, T0, nPoints,forcingdof);

%%%%%%%% Above is forcing setting and set to DynamicalSystem class

clusterRun = true; %% if the script is run on local or cluster.
method="Newmark";
PSDpair=[1,1];

%%
S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
set(S.PSDOptions, 'nPointfilter', 1)
masterModes = [1,2];
S.choose_E(masterModes);
order = 5; % Approximation order
%%
S.ssmSEulerTimeDisp = false; S.tol= 1e-12;
tic
[wss,ssmPSD]=S.extract_PSD(PSDpair, order,'direct heun',[0 25],false);
time_ssm=toc;
disp([num2str(time_ssm),' amount of time in SSM computation'])
%%
char=[num2str(nDiscretization),'QCmodel',num2str(epsilon),'T',num2str(T0),'nP',num2str(log(nPoints)/log(2)),'.mat'];
save(char)
end
