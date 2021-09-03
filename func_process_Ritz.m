function [r_covered, id_cluster_Ritz, n_Ritz_max_target_update] = func_process_Ritz(Ritz, n_converged_Ritz, n_Ritz_max_target, dim_subspace_kept)
n_Ritz = size(Ritz,1);
id_cluster_Ritz = zeros(n_Ritz, 1);
tol_r_cover = 0;
if n_converged_Ritz == 0 % no Ritz value converged, reorder Ritz value by their amplitude
    [~, id_cluster_Ritz] = sort(abs(Ritz), 'descend');
    r_covered = 0;
    n_Ritz_max_target_update = n_Ritz_max_target;
else % Ritz value converged
    Ritz_converge = Ritz(1:n_converged_Ritz);
    min_Ritz_converge = min(abs(Ritz_converge));
    Ritz_unconverge_valid = Ritz((n_converged_Ritz + 1) : dim_subspace_kept);
    idx_Ritz_delay = abs(Ritz_unconverge_valid) > min_Ritz_converge;
    if any(idx_Ritz_delay) % evaluate coverage of converged Ritz values
        a_Ritz_cmp = max(abs(Ritz_unconverge_valid(idx_Ritz_delay)));
        idx_Ritz_converge_sure = abs(Ritz_converge) > a_Ritz_cmp;
        if any(idx_Ritz_converge_sure)
            r_covered = max([1./min(abs(Ritz_converge(idx_Ritz_converge_sure))), tol_r_cover]);
        else
            r_covered = 0;
        end
    else
        r_covered = max([1./min(abs(Ritz_converge)),tol_r_cover]);
    end
    id_cluster_Ritz(1:n_converged_Ritz) = 3;
    id_cluster_Ritz(n_converged_Ritz + find(idx_Ritz_delay)) = 2;
    id_cluster_Ritz(n_converged_Ritz + find(~idx_Ritz_delay)) = 1;
    
    n_Ritz_max_target_update = max([n_converged_Ritz + sum(idx_Ritz_delay), n_Ritz_max_target]); % updation of subspace dimension is an important technique
end
end