function [Fstochastic] = generate_stochastic(obj)
%%% description: generate one particular realization, time history...
    % This function generates one particular sequence of random process on a given time
    % span. There are two ways, one is to generate the realization directly
    % from desired PSD function as an input, while the other one is to
    % indirectly generate PSD from a linear system subjected to a white
    % noise. The latter one has the idea from that the linear system with
    % white noise input can be identified analytically. From this way, we
    % could generate the desired PSD.
    
    % Both methods are done in the second order form
    
    % When the method is direct generation, the input should include a
    % vector of \omega as well as a vector of PSD magnitude values, those
    % are deposited in the object's property 'samplePSD'. When the method
    % is indirect, the input should include linear system's parameters: M C
    % K matrices as well as the white noise magnitude.

    %%%%% note that forcing degree of freedom is essential for finite
    %%%%% element method
    
dim=obj.n;
T0=obj.timeSpan;
N=obj.nPoints;
forcingdof=obj.forcingdof;
    switch obj.SSOptions.ssMethod
        case 'indirect'
            PSD = obj.filterPSD;
            M = PSD.Mz;
            C = PSD.Cz;
            K = PSD.Kz;
            WhiteNoise_S = PSD.S;
            T = linspace(0,T0,N+1); obj.timevector = T;
            m = length(M);
            % SDE method
                Zx=zeros(m,N+1);
                Zv=zeros(m,N+1);
                % sub functions generating forcing realizations/ also forcing
                % directly from PSD
                sigma=sqrt(WhiteNoise_S*2*pi); % variance of white noise with intensity S
                % right now we only deal with static noise
                dW = zeros(m,1);   
                detT=T0/N;
            for i=1:N
                detu = sigma*randn*sqrt(detT);
                dW(end)=detu;
                
                
                Zxhat = (M+1/2*C*detT+1/4*K*detT^2)\((M+1/2*C*detT)*...
                    Zx(:,i)+1/2*detT*(M*Zv(:,i)+1/2*dW));
                Zvhat=(M+1/2*C*detT)\(M*Zv(:,i)-1/2*detT*K*Zxhat+1/2*dW);
                
                % update
                Zx(:,i+1)=Zx(:,i)+detT*Zvhat;
                Zv(:,i+1)=Zv(:,i)-M\(C*Zvhat*detT+...
                    K*Zxhat*detT-dW);
                
            end
                % assigning stochastic force
                Fext=zeros(dim,N+1);
                for j=1:length(forcingdof)
                    l=forcingdof(j);  
                    Fext(l,:)=Zv(1,:); % taking velocity
                end
                %%%
        case 'direct' %% need Phi and omega vectors
            samplePSD=obj.samplePSD;
            Phi=samplePSD(1,:);
            omega=samplePSD(2,:);
            %%%% numebr of sampling points 2^N assumed to be even
            N_sample = length(Phi);
            if (-1)^N_sample<0
                N_sample = N_sample -1;
            end
            %%%% omega assumed to be equally spaced sampling
            w0=omega(2)-omega(1);
            %%%% sampling frequency, time spacing
            T_span=2*pi/w0; T=linspace(0,T_span,N_sample+1);
            obj.sde.timevector = T;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% NOTE: only half of sampled points are used
            Cx=zeros(1,N_sample/2+1);
            for k=1:N_sample/2
                Cx(k)=sqrt(Phi(k)*w0); % amplitude
                theta=2*pi*rand;    %%generate uniformly distributed phase
                Cx(k)=Cx(k)*exp(1i*theta);
            end

            Cxx=zeros(1,N_sample+1);

            for i=2:N_sample/2+1
                Cxx(i)=Cx(i);
            end
            for i=0:N_sample/2-1
                Cxx(N_sample/2+2+i)=conj(Cx(N_sample/2+1-i));
            end
            x=ifft(Cxx)*length(Cxx);
            % assigning stochastic force
                Fext=zeros(dim,length(x));
                for j=1:length(forcingdof)
                    l=forcingdof(j);  
                    Fext(l,:)=x;
                end
                %%%
    end
    
                %%% deleting first zero elements
                for k=1:N+1
                    for f=1:dim
                        if Fext(f,k)~=0
                           stp_i=k;
                            break;
                        end
                    end 
                    
                        if Fext(f,k)~=0
                           stp_i=k;
                            break;
                        end
                end
    Fstochastic=[Fext(:,stp_i:end),zeros(dim,stp_i-1)];
end
