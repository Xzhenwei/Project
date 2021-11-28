function [w,PSD_a]=monte_carlo_average(obj,method,PSDpair,nRealization)

    if nRealization>1
        euler = parcluster('local');
        pool = parpool(euler);

        X = 0; w = 0;
        parfor i=1:nRealization
            [ws,outputPSD] = sde_solver(obj,method,PSDpair);
            X = X+outputPSD;
            w = w + ws;
            disp(['number of realizations left:', num2str(nRealization-i)])
        end


        PSD_a = X/nRealization;
        w = w/nRealization;

        pool.delete()
    else

        [w,PSD_a]=obj.sde_solver(method,PSDpair);


    end

end

