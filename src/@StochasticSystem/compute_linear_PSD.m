function [w,linear_analytic]=compute_linear_PSD(obj,PSDpair,freq_range)
% linear_analytic. This function computes the correspoding linear response
% PSD with the input PSD along with its frequency domain vector omega.

M = obj.M; C = obj.C; K = obj.K;
forcingdof = obj.forcingdof;

N = obj.nPoints;
T0 = obj.timeSpan * 1;
w = (1:N+1)*2*pi/T0;
nOutput = size(PSDpair,1);
linear_analytic = zeros(nOutput , N+1);

if obj.n > 30

        euler = parcluster('local');
        pool = parpool(euler);

    for i = 1:nOutput
    pp1 = PSDpair(i,1);
    pp2 = PSDpair(i,2);

        switch obj.SSOptions.ssMethod
            case 'indirect'
                PSD = obj.filterPSD;
                Mz = PSD.Mz;
                Cz = PSD.Cz;
                Kz = PSD.Kz;
                S = PSD.S;
                G = PSD.G;
        %%%% use par loop
        %%%% check the inverse
                parfor j = 1:N + 1
        %             Hw_z = inv(-w(j)^2*Mz+1i*w(j)*Cz+Kz)*w(j)*1i; %% taking the velocity of the auxilliary system
                    if w(j) < 2*max(freq_range)
                        Hw_z = inv(-w(j)^2*Mz+1i*w(j)*Cz+Kz); %% taking the displacement of the auxilliary system
%                         Hw_z = inv(-w(j)^2*Mz+1i*w(j)*Cz+Kz)*w(j)*1i; %% taking the velocity of the auxilliary system
                        Phi_F = G*G'*S;
                        Zj = Hw_z*Phi_F*(Hw_z');

                        Hw = (-w(j)^2*M+1i*w(j)*C+K);
                        Z_full = Hw\Zj/(Hw'); 

                        linear_analytic(i,j) = norm(Z_full(pp1,pp2));
                    else 
                        linear_analytic(i,j) = 0;
                    end
                end
            case 'direct'
                w = obj.samplePSD(2,:);
                nOmega = length(w);
                Zj = zeros(obj.n,obj.n);
                for j = 1:nOmega 
                    for l = 1:length(forcingdof)
                        Zj(forcingdof(l),forcingdof(l)) = obj.samplePSD(1,j);
                    end
                    if w(j) < 2*max(freq_range)
                        Hw = (-w(j)^2*M+1i*w(j)*C+K);
                        Z_full = Hw\Zj/(Hw'); 

                        linear_analytic(i,j) = norm(Z_full(pp1,pp2));
                    else 
                        linear_analytic(i,j) = 0;
                    end
                end
        
            otherwise
                error('please specify an input PSD')

        end
        pool.delete()
    end
%%%%%%%%%%%%%%%
else
    for i=1:nOutput
        pp1 = PSDpair(i,1);
        pp2 = PSDpair(i,2);

        switch obj.SSOptions.ssMethod
            case 'indirect'
                PSD = obj.filterPSD;
                Mz = PSD.Mz;
                Cz = PSD.Cz;
                Kz = PSD.Kz;
                S = PSD.S;
                G = PSD.G;
        %%%% use par loop
        %%%% check the inverse
                for j = 1:N + 1
        %             Hw_z = inv(-w(j)^2*Mz+1i*w(j)*Cz+Kz)*w(j)*1i; %% taking the velocity of the auxilliary system
                    if w(j) < 2*max(freq_range)
                        Hw_z = inv(-w(j)^2*Mz+1i*w(j)*Cz+Kz); %% taking the displacement of the auxilliary system
%                         Hw_z = inv(-w(j)^2*Mz+1i*w(j)*Cz+Kz)*w(j)*1i; %% taking the velocity of the auxilliary system
                        Phi_F = G*G'*S;
                        Zj = Hw_z*Phi_F*(Hw_z');

                        Hw = (-w(j)^2*M+1i*w(j)*C+K);
                        Z_full = Hw\Zj/(Hw'); 

                        linear_analytic(i,j) = norm(Z_full(pp1,pp2));
                    else 
                        linear_analytic(i,j) = 0;
                    end
%                     disp([num2str((N-j-1)/N) , ' %'])
                end
            case 'direct'
                w = obj.samplePSD(2,:);
                nOmega = length(w);
                Zj = zeros(obj.n,obj.n);
                for j = 1:nOmega 
                    for l = 1:length(forcingdof)
                        Zj(forcingdof(l),forcingdof(l)) = obj.samplePSD(1,j);
                    end
                    if w(j) < 2*max(freq_range)
                        Hw = (-w(j)^2*M+1i*w(j)*C+K);
                        Z_full = Hw\Zj/(Hw'); 

                        linear_analytic(i,j) = norm(Z_full(pp1,pp2));
                    else 
                        linear_analytic(i,j) = 0;
                    end
                end

            otherwise
                error('please specify an input PSD')

        end

    end
%%%%%%%%%%%%%%%
end
xlab = '$\omega$ Frequency';
ylab = 'Power';
linewidth = 1.8;

    figure
    x1 = PSDpair(i,1);
    x2 = PSDpair(i,2);
            
    plot(w,linear_analytic(i,:),'linewidth',linewidth,'DisplayName',...
                'Linear Analytic PSD');
    hold on
    title (['PSD of Dof ',num2str(x1),' and ',num2str(x2)])
    add_labels(xlab,ylab)
    lgd = legend();
    set(lgd,'Interpreter','latex','Location','best');
    grid on
    xlim([min(freq_range) max(freq_range)])
    hold off

end


function add_labels(xlab,ylab)
xlabel(xlab,'Interpreter','latex');
ylabel(ylab,'Interpreter','latex');
set(gca,'FontSize',14);
grid on, axis tight; legend boxoff;
end
