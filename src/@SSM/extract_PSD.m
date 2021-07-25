function [w,Gs] = extract_PSD(obj, PSDpair, ORDER, method,clusterRun)

% obj.System.input_PSD();

nRealization = obj.System.nRealization;

Gs=zeros(numel(ORDER),obj.System.nPoints+1);

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

    [w, X_l] = compute_analyticPSD(obj,PSDpair);
    
    num_points = obj.System.nPoints;
    
    nOutput=size(PSDpair,1);
    
    for k=1:nOutput
        Gs = zeros(num_points+1,1)'; 

        for i=3:num_points+1
            Gs(k,i)= interp1(wss,Gzz(k,:),w(i)) +X_l(k,i);
        end
    
    end


end
% figrue
% plot(w,Gs,'linewidth',1.5)
% xlim([0 10])


end