function [w,PSD_a]=monte_carlo_average(obj,method,PSDpair,nRealization,clusterRun)

if nRealization>1

    if clusterRun
        euler = parcluster('local');
        pool = parpool(euler,24);
    else
        pool = parpool('local',2);
    end

    [w,X] = obj.sde_solver(method,PSDpair);

    parfor i=1:nRealization-1
        [~,outputPSD] = obj.sde_solver(method,PSDpair);
        X = X+outputPSD;
    end
    
    
    
    PSD_a = X/nRealization;


    pool.delete()
else
    
    [w,PSD_a]=obj.sde_solver(method,PSDpair);
    
    
end



end

