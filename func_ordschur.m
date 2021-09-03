function [ut, ht] = func_ordschur(h, n_eig_converged, dim_subspace_kept, r_circle)
    % order rule:
    % part 1: converged Ritz | part2 : un-converged Ritz within r_circle & within kept subspace |
    % part 3: un_converged Ritz outside r_circle & within kept subspace |
    % part 4: Ritz outside kept subspace
    n = size(h,1);
    d = 1./abs(diag(h));
    
    id_cluster = zeros(n,1);
    
    idx_Ritz_converged = ((1:n)' <= n_eig_converged);
    idx_Ritz_kept = ((1:n)' <= dim_subspace_kept);
    
    idx_Ritz_reorder_in_circle = (idx_Ritz_kept & (~idx_Ritz_converged) & d <= r_circle);
    n_Ritz_reorder_in_circle = sum(idx_Ritz_reorder_in_circle);
    
    idx_Ritz_reorder_out_circle = (idx_Ritz_kept & (~idx_Ritz_converged) & d > r_circle);
    n_Ritz_reorder_out_circle = sum(idx_Ritz_reorder_out_circle);
    
    if n_Ritz_reorder_out_circle > 0
        id_cluster(idx_Ritz_reorder_out_circle) = 1;
    end    
    
    if n_Ritz_reorder_in_circle > 0
%         [~, idx_sort_diag] = sort(d(idx_Ritz_reorder_in_circle), 'descend');
        id_cluster(idx_Ritz_reorder_in_circle) = 2;
    end
    id_cluster(idx_Ritz_converged) = 3;
    
    [ut, ht] = ordschur(eye(size(h,1)), h, id_cluster);   
end