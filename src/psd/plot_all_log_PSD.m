function plot_all_log_PSD(omega,Gxx,ORDER,PSDpair,freq_range, islog)
colors = get(0,'defaultaxescolororder');
w = omega.w;
w_linear = omega.linear;
Gss = Gxx.Gss;
linear_analytic = Gxx.linear_analytic;

G_full = Gxx.full;
linear_full = Gxx.linear_full;


xlab = '$\omega$ Frequency';
ylab = 'Power(dB)';
linewidth = 1.8;
nOutDof = size(Gss,1)/numel(ORDER);

    for i=1:nOutDof
        figure
        x1 = PSDpair(i,1);
        x2 = PSDpair(i,2);
        if islog
            
%             for j = 1:numel(ORDER)
%                 order = ORDER(j);
% 
%                 semilogy(w,Gss(nOutDof*(j-1)+i,:),'Color',colors(j,:),'linewidth',linewidth,'DisplayName',...
%                     strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ PSD'));
%                 hold on
%             end
% 
%             semilogy(w_linear,linear_analytic(i,:),'Color',colors(j+1,:),'linewidth',linewidth,...
%                 'DisplayName','linear system analytic')
% 
%             semilogy(w,linear_full(i,:),'Color',colors(j+2,:),'linewidth',linewidth,...
%                 'DisplayName','Full linear system response')
%             semilogy(w,G_full(i,:),'Color',colors(j+3,:),'linewidth',linewidth,...
%                 'DisplayName','Full nonlinear system response')
%             add_labels(xlab,'Power')
        
            for j = 1:numel(ORDER)
                order = ORDER(j);

                plot(w,20*log(Gss(nOutDof*(j-1)+i,:)),'Color',colors(j,:),'linewidth',linewidth,'DisplayName',...
                    strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ PSD'));
                hold on
            end

            plot(w_linear,20*log(linear_analytic(i,:)),'Color',colors(j+1,:),'linewidth',linewidth,...
                'DisplayName','linear system analytic')

            plot(w,20*log(linear_full(i,:)),'Color',colors(j+2,:),'linewidth',linewidth,...
                'DisplayName','Full linear system response')
            plot(w,20*log(G_full(i,:)),'Color',colors(j+3,:),'linewidth',linewidth,...
                'DisplayName','Full nonlinear system response')
            add_labels(xlab,ylab)
        else
            for j = 1:numel(ORDER)
                order = ORDER(j);

                plot(w,Gss(nOutDof*(j-1)+i,:),'Color',colors(j,:),'linewidth',linewidth,'DisplayName',...
                    strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ PSD'));
                hold on
            end

            plot(w_linear,linear_analytic(i,:),'Color',colors(j+1,:),'linewidth',linewidth,...
                'DisplayName','linear system analytic')

            plot(w,linear_full(i,:),'Color',colors(j+2,:),'linewidth',linewidth,...
                'DisplayName','Full linear system response')
            plot(w,G_full(i,:),'Color',colors(j+3,:),'linewidth',linewidth,...
                'DisplayName','Full nonlinear system response')
            add_labels(xlab,'Power')
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