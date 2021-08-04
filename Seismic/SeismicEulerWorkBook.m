function SeismicEulerWorkBook(epsilon)
n = 10;
m = 7;
k = 4555;
c = 90;
kappa2 = 0;
kappa3 = 2;

[M,C,K,fnl] = build_model(n,m,c,k,kappa2,kappa3);
                                                    
nRealization=1;
resol = 1;
T0=resol*100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=resol*2^14; %% control the accuracy of numerical differential equation
[filterPSD,forcingdof,stochastic_f] = build_stochasticF(m,epsilon,n);

SS = StochasticSystem();

set(SS,'filterPSD',filterPSD,'linear',true)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl,'gFactor',-m);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS.SSOptions,'ssMethod','indirect')
SS.add_random_forcing(nRealization, T0, nPoints,forcingdof);

clusterRun=true; %% if the script is run on local or cluster.
method="filter ImplicitMidPoint";
PSDpair=[n,n];

firts_res=abs(imag(D(1)));
SS.sdeImpTimeDisp = false;
[w,outputPSD] = SS.monte_carlo_average(method,PSDpair,nRealization,clusterRun);

[w_linear, linear_analytic]=SS.compute_linear_PSD(PSDpair);

S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
% set(S.Options, 'reltol', 0.1,'notation','tensor')
masterModes = [1,2]; 
S.choose_E(masterModes);

order = 5;
freq_range=[0 7]; % depend on res
S.ssmSEulerTimeDisp = false;
[wss,ssmPSD]=S.extract_PSD(PSDpair, order,'filter heun',freq_range,clusterRun);


figure
plot(wss,ssmPSD,'linewidth',1)
hold on
plot(w,outputPSD(1,:),'linewidth',1)
hold on
plot(w_linear,linear_analytic(1,:),'linewidth',1)
xline(firts_res,'-',{'First Resonance'},'linewidth',1.5);
legend('SSM','nonlinear simulation','linear analytic')
xlim([0,7]);
grid on

