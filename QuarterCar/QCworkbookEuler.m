function QCworkbookEuler(epsilon)
%% kappa=1;
n=2; kappa=2e4;
% kappa=0;
[M,C,K,fnl,f_0] = build_model(kappa,n);
% LINEAR QUARTER CAR with Dimension 2;
nRealization=48;
T0=50; %% PSD frequency domain resolution is ~ w0= 2*pi* 1/T0
nPoints=2^15; %% control the accuracy of numerical differential equation
[samplePSD,forcingdof,IC,stochastic_f] = build_stochasticF(n,nPoints,epsilon);

SS = StochasticSystem();

set(SS,'samplePSD',samplePSD,'linear',false)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS.SSOptions,'ssMethod','direct','tol',1e-12)
SS.add_random_forcing(nRealization, T0, nPoints,forcingdof);

%%%%%%%% Above is forcing setting and set to DynamicalSystem class

clusterRun=false; %% if the script is run on local or cluster.
method="Newmark";
PSDpair=[1,1];
[V,D,W] = SS.linear_spectral_analysis();
firts_res = abs(imag(D(1)));
%%
SS.sdeImpTimeDisp = true;
tic
[w,outputPSD] = SS.monte_carlo_average(method,PSDpair,nRealization,clusterRun);
time_sde = toc;
disp(['Total number of ',num2str(nRealization),'# realization takes ',...
    num2str(time_sde),' amount of time'])
%%
[w_l,linear_analytic]=SS.compute_linear_PSD(PSDpair,[0,25]);
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
[wss,ssmPSD]=S.extract_PSD(PSDpair, order,'direct heun',[0 25],clusterRun);
time_ssm=toc;
disp([num2str(time_ssm),' amount of time in SSM computation'])
%%
tic
[w_galerkin, PSD_galerkin] =  SS.monte_carlo_galerkin(method, PSDpair);
time_galerkin = toc;
disp(['Total number of ',num2str(nRealization),'# realization in Galerkin projection takes ',...
    num2str(time_galerkin),' amount of time'])
%%
char=['QCmodel',num2str(epsilon),'T',num2str(T0),'nP',num2str(log(nPoints)/log(2)),'.mat'];
save(char)
end
