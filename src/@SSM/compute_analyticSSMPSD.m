function [w, X_l] = compute_analyticSSMPSD(obj,PSDpair,freq_range ,clusterRun)
N = obj.System.nPoints;
T0 = obj.System.timeSpan;
w = (1:N+1)*1/T0*2*pi;
Wnode = obj.E.adjointBasis';
    B=obj.System.B; M=obj.System.M; C=obj.System.C; K=obj.System.K; 
    n=obj.System.n;
    G_tf=eye(2*n)-B*obj.E.basis*Wnode;
    G11=G_tf(1:n,1:n);
    X_l=zeros(size(PSDpair,1),N+1);
    
    if clusterRun
        euler = parcluster('local');
        pool = parpool(euler,24);
    else
        pool = parpool('local',2);
    end
    
for i=1:size(PSDpair,1)
    
switch obj.System.SSOptions.ssMethod
    case 'indirect'
        PSD = obj.System.filterPSD;
        Mz = PSD.Mz;
        Cz = PSD.Cz;
        Kz = PSD.Kz;
        S = PSD.S;
        G = PSD.G;

        parfor j = 1:N + 1
            if w(j) < 2*max(freq_range)
                Hw_z = (-w(j)^2*Mz+1i*w(j)*Cz+Kz); %% taking the 
                %%% displacement of the auxilliary system
                Phi_F = G*G'*S;
                Zj = Hw_z\Phi_F/(Hw_z');

                Hw = (-w(j)^2*M+1i*w(j)*C+K);
                Z_full = Hw\G11*Zj*G11'/(Hw');

                X_l(i,j) = norm(Z_full(PSDpair(i,1),PSDpair(i,2)));
            else
                X_l(i,j) = 0;
            end
%             disp([num2str((1-i/N)*100),' %'])
        end
    case 'direct'
       
%         obj.input.omega=obj.samplePSD(2,:); STILL NEED TO FIX FOR DIRECT
%         obj.input.PSD=obj.samplePSD(1,:);
        
    otherwise
        error('please specify an input PSD')
        
end
end
    pool.delete()
    


end