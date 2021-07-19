function [PSD, stochastic_f] = build_stochasticF(nPoints,epsilon)
Mi=[5,0;0,5];
Ci=[110,10;10,110];
Ki=[20,0;0,20];
Si=epsilon^2*1e3; % white noise intensity
PSD.Mz=Mi;
PSD.Cz=Ci;
PSD.Kz=Ki;
PSD.S=Si; % random forcing parameters subject to intensity
%%%%%%%%
omega=linspace(0,200,nPoints+1);

PSD1 = epsilon^2 * 50 ./ (omega.^5) .* exp(-400./(omega.^4)) ; PSD1(1) = 0;

samplePSD=[PSD1;omega];


stochastic_f=true; % automatically to be true unless specified
