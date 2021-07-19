function [wss,Gss] = compute_ssmPSD(obj, PSDpair, ORDER, method,clusterRun)

obj.System.input_PSD();

nRealization = obj.System.nRealization;
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
        [w,outputPSD] = obj.extract_PSD( PSDpair, ORDER, method);
        Gss = Gss+outputPSD;
        wss=wss+w;
    end
    Gss = Gss/nRealization;
    wss = wss/nRealization;
    pool.delete()
else
    
    [wss,Gss] = obj.extract_PSD( PSDpair, ORDER, method);
end


end