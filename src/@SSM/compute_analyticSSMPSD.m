function [w, X_l] = compute_analyticSSMPSD(obj,PSDpair,freq_range)
N = obj.System.nPoints;
T0 = obj.System.timeSpan;
inputForcingType = obj.System.inputForcingType;
w = (1:N+1)*1/T0*2*pi; 
Wnode = obj.E.adjointBasis';
    B=obj.System.B; M=obj.System.M; C=obj.System.C; K=obj.System.K; 
    n=obj.System.n;
    G_tf=eye(2*n)-B*obj.E.basis*Wnode;
    G11=G_tf(1:n,1:n);
    X_l=zeros(size(PSDpair,1),N+1);
    
    
for i=1:size(PSDpair,1)
    pp1 = PSDpair(i,1);
    pp2 = PSDpair(i,2);
switch obj.System.SSOptions.ssMethod
    case 'indirect'
        PSD = obj.System.filterPSD;
        Mz = PSD.Mz;
        Cz = PSD.Cz;
        Kz = PSD.Kz;
        S = PSD.S;
        G = PSD.G;
        if obj.System.n > 200
            euler = parcluster('local');
            pool = parpool(euler);
        
            parfor j = 1:N + 1
                Zj = 0;
                if w(j) < 10*max(freq_range)
                    Hw_z = (-w(j)^2*Mz+1i*w(j)*Cz+Kz); %% taking the 

                    Phi_F = G*G'*S;
                    switch inputForcingType
                        case 'disp'
                            Zj = Hw_z\Phi_F/(Hw_z'); %%% displacement 
                        case 'vel'
                            Zj = 1i*w(j)*Hw_z\Phi_F/(Hw_z')*1i*w(j);%%% vel
                    end

                    Hw = (-w(j)^2*M+1i*w(j)*C+K);
                    Z_full = Hw\G11*Zj*G11'/(Hw');

                    X_l(i,j) = norm(Z_full(pp1,pp2));
                else
                    X_l(i,j) = 0;
                end
            end
            pool.delete()
        else
            for j = 1:N + 1
                if w(j) < 2*max(freq_range)
                    Hw_z = (-w(j)^2*Mz+1i*w(j)*Cz+Kz); %% taking the 

                    Phi_F = G*G'*S;
                   switch inputForcingType
                        case 'disp'
                            Zjj = Hw_z\Phi_F/(Hw_z'); %%% displacement 
                        case 'vel'
                            Zjj = 1i*w(j)*Hw_z\Phi_F/(Hw_z')*1i*w(j);%%% vel
                   end
                    
                    Hw = (-w(j)^2*M+1i*w(j)*C+K);
                    Z_full = Hw\G11*Zjj*G11'/(Hw');

                    X_l(i,j) = norm(Z_full(pp1,pp2));
                else
                    X_l(i,j) = 0;
                end
%             disp([num2str((1-i/N)*100),' %'])
            end
        end
            
    case 'direct'
        forcingdof = obj.System.forcingdof ;
                w = obj.System.samplePSD(2,:);
                nOmega = length(w);
                Zj = zeros(obj.System.n,obj.System.n);
                for j = 1:nOmega 
                    for l = 1:length(forcingdof)
                        Zj(forcingdof(l),forcingdof(l)) = obj.System.samplePSD(1,j);
                    end
                    if w(j) < 2*max(freq_range)

                        Hw = (-w(j)^2*M+1i*w(j)*C+K);
                        Z_full = Hw\G11*Zj*G11'/(Hw');

                        X_l(i,j) = norm(Z_full(pp1,pp2));
                    else
                        X_l(i,j) = 0;
                    end
                    
                end
        
    otherwise
        error('please specify an input PSD')
        
end
end



end