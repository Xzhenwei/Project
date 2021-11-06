function [x, T] = generate_X_0 (omega,Phi)
% samplePSD=obj.samplePSD;
%             Phi=samplePSD(1,:);
%             omega=samplePSD(2,:);
            %%%% numebr of sampling points 2^N assumed to be even
            N_sample = length(Phi);
            if (-1)^N_sample<0
                N_sample = N_sample -1;
            end
            %%%% omega assumed to be equally spaced sampling
            w0=omega(2)-omega(1);
            %%%% sampling frequency, time spacing
            T_span=2*pi/w0; T=linspace(0,T_span,N_sample+1);

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
end