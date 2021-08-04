function [w,Gss] = extract_PSD(obj, PSDpair, ORDER, method,freq_range, clusterRun)


num_points = obj.System.nPoints;

nRealization = obj.System.nRealization;
nOutput      = size(PSDpair,1);

Gs  = zeros(nOutput,obj.System.nPoints+1);
Gss = zeros(nOutput*numel(ORDER),obj.System.nPoints+1);

colors = get(0,'defaultaxescolororder');
    for j = 1:numel(ORDER)
        order = ORDER(j);
        [W0, R0] = obj.compute_whisker(order);

        disp(['Compute the PSD through SSM of order ',num2str(order)]);
        if nRealization>1

            if clusterRun
                euler = parcluster('local');
                pool = parpool(euler,24);
            else
                pool = parpool('local',2);
            end

        %     [w,Gss] = obj.extract_PSD( PSDpair, ORDER, method);
            Gzz=0; wss=0;
            parfor i=1:nRealization
                [w,outputPSD] = compute_ssmPSD(obj, PSDpair, W0, R0, method);
                Gzz = Gzz+outputPSD;
                wss = wss+w;
            end
            Gzz = Gzz/nRealization;
            wss = wss/nRealization;
            pool.delete()
        else

            [wss,Gzz] = compute_ssmPSD(obj, PSDpair, W0, R0, method);
        end

        [w, X_l] = compute_analyticSSMPSD(obj,PSDpair);




        for k=1:nOutput


            for i=1:num_points+1
                Gs(k,i)= interp1(wss,Gzz(k,:),w(i)) +X_l(k,i);
            end

        end

    Gss(nOutput*(j-1)+1:nOutput*j,:) = Gs;
    
    
    end
    
    
plot_PSD(w,Gss,ORDER,PSDpair,colors,freq_range)


end