function [w,outputPSD]=sde_solver(obj,SDEmethod,PowerSpectralPair) 

N = obj.nPoints;
T0 = obj.timeSpan;
dim = obj.n;

% initial conditions are set to be 0
X0 = zeros(dim,1); V0 = X0;

    if ~isempty(obj.filterPSD)
        PSD = obj.filterPSD;
        Mz = PSD.Mz;
        Cz = PSD.Cz;
        Kz = PSD.Kz;
        WhiteNoise_S = PSD.S;
        StochasticRealization = "indirect";

    elseif ~isempty(obj.samplePSD)
       
       StochasticRealization= "direct";
       
    else
        error('please specify an input PSD')
        
    end



%%
%%% initialization for output PSD
% SpecDensPair format ie: [1,1;2,2;3,3]
nPSDpairs = size(PowerSpectralPair,1);
Cy = zeros(nPSDpairs,N+1);   
        T=linspace(0,T0,N+1); detT=T0/N;
        % SDE
            switch SDEmethod
                case "filter ImplicitMidPoint"
                    G=zeros(obj.n,length(Mz)); %G(end,:)=ones(1,length(Mz));
                    PSD.G=G;
                    [X,V]=implicit_Mid_Point(obj,N,T0,PSD);
                case "filterHeun"
            %Forward Heun's method
                    obj.Fsto = generate_stochastic(obj,StochasticRealization);

                    G=zeros(obj.n,length(Mz)); %G(end,:)=ones(1,length(Mz));
                    PSD.G=G;
                    [X,V]=forward_Heun(obj,N,T0,PSD);
            % Newmark method
                case "Newmark"
                    obj.Fsto = generate_stochastic(obj,StochasticRealization);
                    
                    A0=zeros(dim,1);
                    TI_sto = ImplicitNewmark('timestep',detT,'alpha',0.005,'linear',obj.linear);
                    % Modal linear Residual evaluation function handle
                    Residual_sto = @(q,qd,qdd,t)residual(obj,q,qd,qdd,t);
                    % time integration
                    TI_sto.Integrate(X0,V0,A0,T0,Residual_sto);
                    X=TI_sto.Solution.q;
             end
% SpecDensPair indicates the power spectral density. format ie:
% [1,1;2,2;3,3] meaning we want to calculate PSD of x1, x2,x3 respectively
        
        for j=1:nPSDpairs
            [w,Cy(j,:)]=crossPSDestimator(X(PowerSpectralPair(j,1),:),X(PowerSpectralPair(j,2),:),T);
        end
        
        outputPSD=Cy;
    
    if StochasticRealization== "indirect"
        Z11 = zeros(1, N+1);
        Phi_F = zeros(length(Mz), length(Mz));
        for j = 1:N + 1
         Hw = inv(-w(j)^2*Mz+1i*w(j)*Cz+Kz);
         Phi_F(end,end) = WhiteNoise_S;
         Zj=(-w(j)^2*Mz+1i*w(j)*Cz+Kz)\Phi_F*Hw.';
         Z11(j)=norm(Zj(1,1)); %% this is the input analytical signal/process
        end
        obj.input.omega=w;
        obj.input.PSD=Z11; %% we store the input signal as the last row in outputPSD
    else
        obj.input.omega=obj.samplePSD(2,:);
        obj.input.PSD=obj.samplePSD(1,:);
    end
end
function [w,Gxy]=crossPSDestimator(x,y,t)
%%% Reference: Preunont Andre, Random vibration and Spectral Analysis Chp 12.9
%%%% numebr of sampling points 
N = length(x); 
%%%% time period of generation
T0=max(t);
%%%% sampling frequency, time spacing
f0=1/T0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cx=fft(x)/N;Cy=fft(y)/N;Cys=conj(Cy);

Gxy=T0*(Cx.*Cys)/2/pi;
w=(1:length(Cx))*f0*2*pi;

end