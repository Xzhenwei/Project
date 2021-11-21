function SeismicEulerWorkBook(epsilon)
n = 10;
m = 7;
k = 4555;
c = 90;
kappa2 = 0;
kappa3 = 1.5;

[M,C,K,fnl] = build_model(n,m,c,k,kappa2,kappa3);
nRealization=24;
T0=100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=2^14; %% control the accuracy of numerical differential equation
 %% forcing magnitude
[filterPSD,forcingdof,stochastic_f] = build_stochasticF(m,epsilon,n);
%%%%%%%%
SS = StochasticSystem();

set(SS,'filterPSD',filterPSD,'linear',false)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl,'gFactor',-m);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS.SSOptions,'ssMethod','indirect')
SS.add_random_forcing(nRealization, T0, nPoints,forcingdof);

%%%%%%%% Above is forcing setting and set to DynamicalSystem class
clusterRun = true; %% if the script is run on local or cluster.
method="filter ImplicitMidPoint";
PSDpair=[n,n];
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

%% SSM setting
S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
set(S.PSDOptions, 'nPointfilter', 0.5)
masterModes = [1,2];
S.choose_E(masterModes);
order = 5; % SSM approximation order
%%
freq_range=[0 10];
S.ssmSEulerTimeDisp = false;
tic
[wss,ssmPSD]=S.extract_PSD(PSDpair, order,'filter',freq_range,clusterRun);
time_ssm=toc;
disp([num2str(time_ssm),'s amount of time'])
%% linear
[w_linear, linear_analytic]=SS.compute_linear_PSD(PSDpair,freq_range);
%%

[U, OmegaSq, NOT_CONVERGED] = eigs(sparse(K),sparse(M),2,'smallestabs');
% sort OmegaSq and take the largest one
% disp('OmegaSq') to see if the eigenvalues similar
U = U(:,1) ;
% record time and place in table
euler = parcluster('local');
pool = parpool(euler);
PSD_galerkin = 0; w_galerkin = 0;
tic
parfor i=1:nRealization
    
    [w_g,Gz_g] = galerkin_proj(SS, U, filterPSD, method, PSDpair);
    PSD_galerkin = PSD_galerkin+Gz_g;
    w_galerkin = w_galerkin + w_g;
    disp(['number of realizations left:', num2str(nRealization-i)])
    
end
time_galerkin = toc;
PSD_galerkin = PSD_galerkin/nRealization;
w_galerkin = w_galerkin/nRealization;
pool.delete()

%%
char=['SeismicEpsilon',num2str(epsilon),'.mat'];
save(char,'-mat')

end
