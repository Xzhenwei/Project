function plot_log_PSD(omega,Gxx,ORDER,PSDpair,freq_range,islog)
colors = get(0,'defaultaxescolororder');
w = omega.w;
w_linear = omega.linear;
Gss = Gxx.Gss;
linear_analytic = Gxx.linear_analytic;
xlab = '$\omega$ Frequency';
ylab = 'Power';
linewidth = 1.5;
nOutDof = size(Gss,1)/numel(ORDER);

    for i=1:nOutDof
        figure
        x1 = PSDpair(i,1);
        x2 = PSDpair(i,2);
        
        if islog
            for j = 1:numel(ORDER)
                order = ORDER(j);

                plot(w,20*log(Gss(nOutDof*(j-1)+i,:)),'--','Color',colors(j,:),'linewidth',linewidth,'DisplayName',...
                    strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ PSD'));
                hold on
            end

            plot(w_linear,20*log(linear_analytic(i,:)),'Color',colors(j+2,:),'linewidth',linewidth,...
                'DisplayName','linear system response')
            add_labels(xlab,'Power(dB)')
        else
            %%%
            for j = 1:numel(ORDER)
                order = ORDER(j);

                plot(w,(Gss(nOutDof*(j-1)+i,:)),'--','Color',colors(j,:),'linewidth',linewidth,'DisplayName',...
                    strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ PSD'));
                hold on
            end

            plot(w_linear,(linear_analytic(i,:)),'Color',colors(j+2,:),'linewidth',linewidth,...
                'DisplayName','linear system response')
            add_labels(xlab,ylab)
            %%%
        end
        title (['PSD of Dof (',num2str(x1),', ',num2str(x2),')'])
        
        lgd = legend();
        set(lgd,'Interpreter','latex','Location','best');
        grid on
        xlim([min(freq_range) max(freq_range)])
        hold off
    end
    
end

function add_labels(xlab,ylab)
xlabel(xlab,'Interpreter','latex');
ylabel(ylab,'Interpreter','latex');
set(gca,'FontSize',14);
grid on, axis tight; legend boxoff;
end