nElements = 10;
epsilon = 5e-3;

[M,C,K,fnl,outdof,eMass] = build_model(nElements);
f2 = sptensor([n,n,n]); % coeff of quadratic nonlinearity/ 0 
f3 = sptensor([n,n,n,n]); % coeff of cubic nonlinearity / single cubic spring
fnl = {f2,f3};
%%
n = length(M);
nRealization=20;
resol = 1;
T0=resol*100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=resol*8*2^14; %% control the accuracy of numerical differential equation

[filterPSD, stochastic_f] = build_stochasticF(eMass,n,epsilon);
disp(['Number of degrees of freedom = ' num2str(n)])


                                                       
SS = StochasticSystem();

set(SS,'filterPSD',filterPSD,'linear',true)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl,'gFactor',-eMass);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS.SSOptions,'ssMethod','indirect')
SS.add_random_forcing(nRealization, T0, nPoints,outdof);

%%%%%%%% moving linear beam(look for reference)
clusterRun=true; %% if the script is run on local or cluster.
method="Newmark";
PSDpair=[n-1,n-1];

[V, D, W] = SS.linear_spectral_analysis();
firts_res=abs(imag(D(1)));

SS.sdeImpTimeDisp = false;
[w,outputPSD] = SS.monte_carlo_average(method,PSDpair,nRealization,clusterRun);

S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
% set(S.Options, 'reltol', 0.1,'notation','tensor')
masterModes = [1,2]; 
S.choose_E(masterModes);

order = [5];  
freq_range=[4 7];
S.ssmSEulerTimeDisp = false;
[wss,ssmPSD]=S.extract_PSD(PSDpair, order,'filter heun',freq_range,clusterRun);


[w_linear, linear_analytic]=SS.compute_linear_PSD(PSDpair,freq_range);

% char=['vonEpsilon',num2str(epsilon),'.mat'];
% save(char,'-mat')
%%
plot_all_PSD(w,ssmPSD,outputPSD,linear_analytic,order,PSDpair,freq_range,true)
saveas(gcf,'plotall.png')

figure
plot(w,outputPSD(1,:),'linewidth',1.5,... 
            'DisplayName','full system computation')
hold on
plot(w_linear,linear_analytic(1,:),'linewidth',1.5,...
            'DisplayName','linear system response')

xlim([0 10])
saveas(gcf,'plotfullandlin.png')
%%
save('linearwork.mat')