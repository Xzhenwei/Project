function [w_galerkin, PSD_galerkin] = monte_carlo_galerkin(obj, method, PSDpair)
[U, ~, ~] = eigs(sparse(obj.K),sparse(obj.M),2,'smallestabs');

U = U(:,1) ;
filterPSD = obj.filterPSD;
% record time and place in table
PSD_galerkin = 0; w_galerkin = 0;
if obj.nRealization > 1
    euler = parcluster('local');
    pool = parpool(euler);

    parfor i=1:nRealization

        [w_g,Gz_g] = galerkin_proj(obj, U, filterPSD, method, PSDpair);
        PSD_galerkin = PSD_galerkin+Gz_g;
        w_galerkin = w_galerkin + w_g;
        disp(['number of realizations left:', num2str(nRealization-i)])

    end
    PSD_galerkin = PSD_galerkin/nRealization;
    w_galerkin = w_galerkin/nRealization;
    pool.delete()
else
    [w_galerkin, PSD_galerkin] = galerkin_proj(obj, U, filterPSD, method, PSDpair);
end


end