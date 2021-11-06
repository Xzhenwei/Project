%% kappa=1;
n=2; kappa=2e4;
% kappa=0;
[M,C,K,fnl,f_0] = build_model(kappa,n);
% LINEAR QUARTER CAR with Dimension 2;
nRealization=20;
T0=50; %% PSD frequency domain resolution is ~ w0= 2*pi* 1/T0
nPoints=2^15; %% control the accuracy of numerical differential equation
epsilon=0.1; %% forcing magnitude
[samplePSD,forcingdof,IC,stochastic_f] = build_stochasticF(n,nPoints,epsilon);
DS = DynamicalSystem();

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
[w,outputPSD] = SS.monte_carlo_average(method,PSDpair,nRealization,clusterRun);
%% linear
set(SS,'linear',true)
[w_ln,outputPSD_ln] = SS.monte_carlo_average(method,PSDpair,nRealization,clusterRun);
%%
[w_l,linear_analytic]=SS.compute_linear_PSD(PSDpair,[0,25]);
%% this section plots the result of the full system response and linear analytical
figure
plot(w,outputPSD(1,:),'linewidth',1)
hold on
plot(w_ln,outputPSD_ln(1,:),'linewidth',2)
% hold on
xlim([0,25])
% xlabel('\Omega frequency')
% ylabel('Power Density')
% legend(strcat('PSD solved by ', ": ",method), 'linear analytic')
% legend('boxoff')
% title('Power Spectral Density of X1')
% grid on
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
[wss,ssmPSD]=S.extract_PSD(PSDpair, order,' ',[0 25],false);
time_ssm=toc;
disp([num2str(time_ssm),' amount of time'])

%% Plotting 
clear omega
clear Gxx
omega.w = w_galerkin; omega.linear = w_l; omega.wss = wss;
Gxx.Gfull = PSD_galerkin; Gxx.linear_analytic = linear_analytic; Gxx.Gss = ssmPSD; 
plot_log_PSD(omega,Gxx,order,PSDpair,[0 30],false)
%%
[U, OmegaSq, NOT_CONVERGED] = eigs(sparse(K),sparse(M),2,'smallestabs');
U = U(:,1) ;
% record time and place in table
euler = parcluster('local');
pool = parpool(euler);
PSD_galerkin = 0; w_galerkin = 0;
parfor i=1:nRealization
    
    [w_g,Gz_g] = GalerkinProj (U, SS, M, C, K);
    PSD_galerkin = PSD_galerkin+Gz_g;
    w_galerkin = w_galerkin + w_g;
    disp(['number of realizations left:', num2str(nRealization-i)])
    
end
PSD_galerkin = PSD_galerkin/nRealization;
w_galerkin = w_galerkin/nRealization;
pool.delete()
%%
plot(w_galerkin,PSD_galerkin)
xlim([0,30])
