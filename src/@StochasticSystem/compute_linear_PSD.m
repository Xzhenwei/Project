function [w,linear_analytic]=compute_linear_PSD(obj,PSDpair)
% linear_analytic. This function computes the correspoding linear response
% PSD with the input PSD along with its frequency domain vector omega.

n = obj.n;
M = obj.M; C = obj.C; K = obj.K;
forcingdof = obj.forcingdof;

N = obj.nPoints;
T0 = obj.timeSpan * 1;
w = (1:N+1)*1/T0*2*pi;
nOutput = size(PSDpair,1);
linear_analytic=zeros(nOutput , N+1);

for i=1:nOutput
 

switch obj.SSOptions.ssMethod
    case 'indirect'
        PSD = obj.filterPSD;
        Mz = PSD.Mz;
        Cz = PSD.Cz;
        Kz = PSD.Kz;
        S = PSD.S;
        G = PSD.G;

        for j = 1:N + 1
%             Hw_z = inv(-w(j)^2*Mz+1i*w(j)*Cz+Kz)*w(j)*1i; %% taking the velocity of the auxilliary system
            
            Hw_z = inv(-w(j)^2*Mz+1i*w(j)*Cz+Kz)*1i; %% taking the displacement of the auxilliary system
            Phi_F = G*G'*S;
            Zj = Hw_z*Phi_F*Hw_z.';
            
            Hw = inv(-w(j)^2*M+1i*w(j)*C+K);
            Z_full = (-w(j)^2*M+1i*w(j)*C+K)\Zj*Hw';
            
            linear_analytic(i,j) = norm(Z_full(PSDpair(i,1),PSDpair(i,2)));
        end
    case 'direct'
       
%         obj.input.omega=obj.samplePSD(2,:); STILL NEED TO FIX FOR DIRECT
%         obj.input.PSD=obj.samplePSD(1,:);
        
    otherwise
        error('please specify an input PSD')
        
end
end
