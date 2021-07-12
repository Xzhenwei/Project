function [w,Gss] = compute_ssmPSD(obj, PSDpair, ORDER, method,clusterRun)

nRealization = obj.System.nRealization;
if nRealization>1

    if clusterRun
        euler = parcluster('local');
        pool = parpool(euler,24);
    else
        pool = parpool('local',2);
    end

    [w,Gss] = obj.extract_PSD( PSDpair, ORDER, method);

    parfor i=1:nRealization-1
        [~,outputPSD] = obj.extract_PSD( PSDpair, ORDER, method);
        Gss = Gss+outputPSD;
    end
    Gss = Gss/nRealization;

    pool.delete()
else
    
    [w,Gss] = obj.extract_PSD( PSDpair, ORDER, method);
end


end