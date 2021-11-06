function input_PSD(obj)

N = obj.nPoints;
T0 = obj.timeSpan;
w = (1:N+1)*1/T0*2*pi;

switch obj.SSOptions.ssMethod
    case 'indirect'
        PSD = obj.filterPSD;
        Mz = PSD.Mz;
        Cz = PSD.Cz;
        Kz = PSD.Kz;
        S = PSD.S;
        G = PSD.G;
        Z11 = zeros(1, N+1);

        for j = 1:N + 1
         Hw = inv(-w(j)^2*Mz+1i*w(j)*Cz+Kz)*w(j)*1i; %% taking the velocity of the auxilliary system
         Phi_F = G*G'*S;
         Zj=Hw*Phi_F*Hw.';
         Z11(j)=abs(Zj(1,1)); %% this is the input analytical signal/process
        end
        obj.input.omega=w;
        obj.input.PSD=Z11; %% we store the input signal as the last row in outputPSD
    
    case 'direct'
       
        obj.input.omega=obj.samplePSD(2,:);
        obj.input.PSD=obj.samplePSD(1,:);
        
    otherwise
        error('please specify an input PSD')
        
end
end