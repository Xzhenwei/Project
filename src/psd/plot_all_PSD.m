function plot_all_PSD(w,Gss,outputPSD,linear_analytic,ORDER,PSDpair,freq_range,dispLinear)
colors = get(0,'defaultaxescolororder');

xlab = '$\omega$ Frequency';
ylab = 'Power';
linewidth = 1.8;
nOutDof = size(Gss,1)/numel(ORDER);

    for i=1:nOutDof
        figure
        x1 = PSDpair(i,1);
        x2 = PSDpair(i,2);
        for j = 1:numel(ORDER)
            order = ORDER(j);
            
            plot(w,Gss(nOutDof*(j-1)+i,:),'Color',colors(j,:),'linewidth',linewidth,'DisplayName',...
                strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ PSD'));
            hold on
        end
        plot(w,outputPSD(i,:),'Color',colors(j+1,:),'linewidth',linewidth,...
            'DisplayName','full system computation')
        if dispLinear
        plot(w,linear_analytic(i,:),'Color',colors(j+2,:),'linewidth',linewidth,...
            'DisplayName','linear system response')
        end
        title (['PSD of Dof ',num2str(x1),' and ',num2str(x2)])
        add_labels(xlab,ylab)
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