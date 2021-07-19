function [wss,Gs] = compute_ssmPSD(obj, PSDpair, ORDER, method,clusterRun)

obj.System.input_PSD();

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
        Gss=0; wss=0;
        parfor i=1:nRealization
            [w,outputPSD] = extract_PSD(obj, PSDpair, W0, R0, method);
            Gss = Gss+outputPSD;
            wss = wss+w;
        end
        Gss = Gss/nRealization;
        wss = wss/nRealization;
        pool.delete()
    else

        [wss,Gss] = extract_PSD(obj, PSDpair, W0, R0, method);
    end

Gs(j,:)=Gss;

end

end