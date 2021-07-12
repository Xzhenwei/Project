function [w,Gss] = compute_ssmPSD(obj, PSDpair, ORDER, method,clusterRun)

if obj.System.nRealization>1

    if clusterRun
        euler = parcluster('local');
        pool = parpool(euler,24);
    else
        pool = parpool('local',2);
    end

    [w,Gss] = obj.extract_PSD(obj, PSDpair, ORDER, method);

    parfor i=1:nRealization-1
        [~,outputPSD] = obj.extract_PSD(obj, PSDpair, ORDER, method);
        Gss = Gss+outputPSD;
    end
    Gss = Gss/nRealization;

    pool.delete()
else
    
    [w,Gss] = extract_PSD(obj, PSDpair, ORDER, method);
end


end