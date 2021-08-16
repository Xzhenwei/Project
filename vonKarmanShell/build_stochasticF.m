function [PSD, stochastic_f] = build_stochasticF(outdof, n, epsilon)
Mi=5;
Ci=100;
Ki=20;

Si=epsilon^2*1e4; % white noise intensity
PSD.Mz=Mi;
PSD.Cz=Ci;
PSD.Kz=Ki;
PSD.S=Si; % random forcing parameters subject to intensity

G = sparse(n,1); 

for i=1:length(outdof)
    indof = outdof(i);
    G(indof,1) = 1;
end

PSD.G=G;

stochastic_f=true; % automatically to be true unless specified
