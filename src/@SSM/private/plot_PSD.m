function plot_PSD(w,Gss,ORDER,PSDpair,colors,freq_range)

xlab = '$\omega$ Frequency';
ylab = 'Power';
linewidth = 1.8;
nOutDof = size(Gss,1)/numel(ORDER);

    for i=1:nOutDof
        % title (PSD _ 10 , 10)
        figure
        x1 = PSDpair(i,1);
        x2 = PSDpair(i,2);
        for j = 1:numel(ORDER)
            order = ORDER(j);
            
            plot(w,Gss(nOutDof*(j-1)+i,:),'Color', colors(j,:),'linewidth',linewidth,'DisplayName',...
                strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ PSD'));
            hold on
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