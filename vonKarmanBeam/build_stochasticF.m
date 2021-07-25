function [PSD, stochastic_f] = build_stochasticF(m, n, epsilon)
Mi=5;
Ci=100;
Ki=20;

Si=epsilon^2*1e2; % white noise intensity
PSD.Mz=Mi;
PSD.Cz=Ci;
PSD.Kz=Ki;
PSD.S=Si; % random forcing parameters subject to intensity

G=zeros(n,1); G(2:3:end)=1;

PSD.G=-m*G;


stochastic_f=true; % automatically to be true unless specified
