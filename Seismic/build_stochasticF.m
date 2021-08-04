function [PSD,forcingDof,stochastic_f] = build_stochasticF(m,epsilon,n)
Mi=[5,0;0,5];
Ci=[110,10;10,110];
Ki=[20,0;0,20];

Mi=5;
Ci=100;
Ki=20;

Si=epsilon^2*1e6; % white noise intensity
PSD.Mz=Mi;
PSD.Cz=Ci;
PSD.Kz=Ki;
PSD.S=Si; % random forcing parameters subject to intensity
PSD.G=-m*ones(n,1);

forcingDof=(1:n)';

stochastic_f=true; % automatically to be true unless specified
