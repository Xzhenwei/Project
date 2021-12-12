function [samplePSD,forcingDof,stochastic_f] = build_stochasticF(n,nPoints,epsilon,kt,ct)

omega=linspace(0,100,nPoints+1);

dispPSD = epsilon^2 * 50 ./ (omega.^5) .* exp(-400./(omega.^4)) ; dispPSD(1) = 0;
velPSD = dispPSD.*omega.^2;
samplePSD=[kt^2*dispPSD + ct^2*velPSD; omega];

forcingDof=[n];
% PSD.G(forcingDof,1)=1;

stochastic_f=true; % automatically to be true unless specified
