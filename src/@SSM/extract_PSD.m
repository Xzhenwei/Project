function [w,Gss] = extract_PSD(obj, PSDpair, ORDER, method,freq_range)


num_points = obj.System.nPoints;

nRealization = obj.System.nRealization;
nOutput = size(PSDpair,1);

Gs  = zeros(nOutput,obj.System.nPoints+1);
Gss = zeros(nOutput*numel(ORDER),obj.System.nPoints+1);

colors = get(0,'defaultaxescolororder');
    for j = 1:numel(ORDER)
        order = ORDER(j);
        [W0, R0] = obj.compute_whisker(order);

        disp(['Compute the PSD through SSM of order ',num2str(order)]);
        
%         [w_l, X_l] = compute_analyticSSMPSD(obj,PSDpair,freq_range,clusterRun);
%         obj.w_l = w_l;
%         obj.X_l = X_l;
        if nRealization>1

            euler = parcluster('local');
            pool = parpool(euler);

            Gzz = 0; wss = 0;
            parfor i=1:nRealization
                [w,outputPSD] = compute_ssmPSD(obj, PSDpair, W0, R0, method);
                Gzz = Gzz+outputPSD;
                wss = wss+w;
                disp(['number of realizations left:', num2str(nRealization-i)])
            end
            Gzz = Gzz/nRealization;
            wss = wss/nRealization;
            pool.delete()
        else

            [wss,Gzz] = compute_ssmPSD(obj, PSDpair, W0, R0, method);
        end

         
%         X_l = zeros(size(Gzz)); w_l = wss;
         w = (1:num_points+1)*1/obj.System.timeSpan*2*pi; 


        for k = 1:nOutput
            
            Gs(k,:) = interp1(wss,Gzz(k,:),w);

        end

    Gss(nOutput*(j-1)+1:nOutput*j,:) = Gs;
    
    
    end
    
    
plot_PSD(w,Gss,ORDER,PSDpair,colors,freq_range)


end