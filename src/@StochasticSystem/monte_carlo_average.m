function [w,PSD_a]=monte_carlo_average(obj,method,PSDpair,nRealization)
assert(nRealization > 1, ' Please use sde.solver for one simulation')

pool = parpool('local',2);

[w,X] = obj.sde_solver(method,PSDpair);

parfor i=1:nRealization-1
    [~,outputPSD] = obj.sde_solver(method,PSDpair);
    X = X+outputPSD;
end
PSD_a = X/nRealization;


pool.delete()

end

