function [w, X_l] = compute_analyticSSMPSD(obj,PSDpair)
N = obj.System.nPoints;
T0 = obj.System.timeSpan;
w = (1:N+1)*1/T0*2*pi;
Wnode = obj.E.adjointBasis';
    B=obj.System.B; M=obj.System.M; C=obj.System.C; K=obj.System.K; 
    n=obj.System.n;
    G_tf=eye(2*n)-B*obj.E.basis*Wnode;
    G11=G_tf(1:n,1:n);
    X_l=zeros(size(PSDpair,1),N+1);
    
for i=1:size(PSDpair,1)
switch obj.System.SSOptions.ssMethod
    case 'indirect'
        PSD = obj.System.filterPSD;
        Mz = PSD.Mz;
        Cz = PSD.Cz;
        Kz = PSD.Kz;
        S = PSD.S;
        G = PSD.G;

        for j = 1:N + 1
            Hw_z = inv(-w(j)^2*Mz+1i*w(j)*Cz+Kz)*w(j)*1i; %% taking the velocity of the auxilliary system
            Phi_F = G*G'*S;
            Zj = Hw_z*Phi_F*Hw_z.';
            
            Hw = inv(-w(j)^2*M+1i*w(j)*C+K);
            Z_full = (-w(j)^2*M+1i*w(j)*C+K)\G11*Zj*G11'*Hw';
            
            X_l(i,j) = norm(Z_full(PSDpair(i,1),PSDpair(i,2)));
        end
    case 'direct'
       
%         obj.input.omega=obj.samplePSD(2,:); STILL NEED TO FIX FOR DIRECT
%         obj.input.PSD=obj.samplePSD(1,:);
        
    otherwise
        error('please specify an input PSD')
        
end
end
    
    


end