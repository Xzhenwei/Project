%% Shallow-curved shell structure with geometric nonlinearities
% Finite element model used in the following reference:
% 
% Jain, S., & Tiso, P. (2018). Simulation-free hyper-reduction for geometrically 
% nonlinear structural dynamics: a quadratic manifold lifting approach. _Journal 
% of Computational and Nonlinear Dynamics_, _13_(7), 071003. <https://doi.org/10.1115/1.4040021 
% https://doi.org/10.1115/1.4040021>
% 
% Finite element code taken from the following package:
% 
% Jain, S., Marconi, J., Tiso P. (2020). YetAnotherFEcode (Version v1.1). Zenodo. 
% <http://doi.org/10.5281/zenodo.4011282 http://doi.org/10.5281/zenodo.4011282>
% 
clear all; close all; clc
% % run ../../install.m
% % change to example directory
% exampleDir = fileparts(mfilename('fullpath'));
% cd(exampleDir)
%% 
% *system parameters*

nDiscretization = 10; % Discretization parameter 30 (#DOFs is proportional to the square of this number)
epsilon = 50;
%% generate model

[M,C,K,fnl,fext,outdof,out] = build_model(nDiscretization);
n = length(M); % number of degrees of freedom
disp(['Number of degrees of freedom = ' num2str(n)])
disp(['Phase space dimensionality = ' num2str(2*n)])
%%
[filterPSD, stochastic_f] = build_stochasticF(outdof,n,epsilon);
%% Dynamical system setup 
% We consider the forced system
% 
% $$\mathbf{M}\ddot{\mathbf{x}}+\mathbf{C}\dot{\mathbf{x}}+\mathbf{K}\mathbf{x}+\mathbf{f}(\mathbf{x},\dot{\mathbf{x}})=\epsilon\mathbf{f}^{ext}(\mathbf{\Omega}t),$$
% 
% which can be written in the first-order form as 
% 
% $$\mathbf{B}\dot{\mathbf{z}}	=\mathbf{Az}+\mathbf{F}(\mathbf{z})+\epsilon\mathbf{F}^{ext}(\mathbf{\phi}),\\\dot{\mathbf{\phi}}	
% =\mathbf{\Omega}$$
% 
% where
% 
% $\mathbf{z}=\left[\begin{array}{c}\mathbf{x}\\\dot{\mathbf{x}}\end{array}\right],\quad\mathbf{A}=\left[\begin{array}{cc}-\mathbf{K} 
% & \mathbf{0}\\\mathbf{0} & \mathbf{M}\end{array}\right],\mathbf{B}=\left[\begin{array}{cc}\mathbf{C} 
% & \mathbf{M}\\\mathbf{M} & \mathbf{0}\end{array}\right],\quad\quad\mathbf{F}(\mathbf{z})=\left[\begin{array}{c}\mathbf{-\mathbf{f}(\mathbf{x},\dot{\mathbf{x}})}\\\mathbf{0}\end{array}\right],\quad\mathbf{F}^{ext}(\mathbf{z},\mathbf{\phi})=\left[\begin{array}{c}\mathbf{f}^{ext}(\mathbf{\phi})\\\mathbf{0}\end{array}\right]$.

SS = StochasticSystem();   
set(SS,'filterPSD',filterPSD,'linear',false)
set(SS,'M',M,'C',C,'K',K,'fnl',fnl);
set(SS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS.SSOptions,'ssMethod','indirect')
% set(DS.Options,'Emax',5,'Nmax',10,'notation','tensor')
nRealization=20;
T0=100; %% PSD frequency domain resolution is ~ 1/T0
nPoints=2^14; %% control the accuracy of numerical differential equation
%% 
% We assume periodic forcing of the form
% 
% $$\mathbf{f}^{ext}(\phi) = \mathbf{f}_0\cos(\phi)=\frac{\mathbf{f}_0}{2}e^{i\phi} 
% + \frac{\mathbf{f}_0}{2}e^{-i\phi}  $$
% 
%%%%%%%% Above is forcing setting and set to DynamicalSystem class
clusterRun=false; %% if the script is run on local or cluster.
method="filter ImplicitMidPoint";
PSDpair=[out,out];

SS.add_random_forcing(nRealization, T0, nPoints,outdof);

%% Linear Modal analysis and SSM setup

[V,D,W] = SS.linear_spectral_analysis();
firts_res = abs(imag(D(1)));
%% 
% *Choose Master subspace (perform resonance analysis)*

S = SSM(SS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
% set(S.Options, 'reltol', 0.1,'notation','tensor')
masterModes = [1,2]; 
S.choose_E(masterModes);
%% PSD using SSMs
% Obtaining *PSD* in reduced-polar coordinate

order = 5; % Approximation order
%% 
S.ssmSEulerTimeDisp = true;

freq_range=[145 155];
tic
[wss,ssmPSD]=S.extract_PSD(PSDpair, order,'filter',freq_range,clusterRun);
time_ssm=toc;
disp([num2str(time_ssm),' amount of time'])
%% Linear part
SS_l = StochasticSystem();   fnl_l = {sptensor([n,n,n]),sptensor([n,n,n,n])};
set(SS_l,'filterPSD',filterPSD,'linear',true)
set(SS_l,'M',M,'C',C,'K',K,'fnl',fnl_l);
set(SS_l.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(SS_l.SSOptions,'ssMethod','indirect')
SS_l.add_random_forcing(nRealization, T0, nPoints,outdof);
S_l = SSM(SS_l);
set(S_l.Options, 'reltol', 0.1,'notation','multiindex')
S_l.choose_E(masterModes);
S_l.ssmSEulerTimeDisp = false;
[w_l,ssmPSD_l]=S_l.extract_PSD(PSDpair, order,'filter heun',freq_range,clusterRun);
%%
% [w_linear, linear_analytic]=SS.compute_linear_PSD(PSDpair,freq_range);
% 
plot_ssm_lin_PSD(log(wss),20*log(ssmPSD),20*log(ssmPSD_l),order,PSDpair,log(freq_range),true)

