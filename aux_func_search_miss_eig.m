function [idx_eig_target_loss, max_mismatch] = aux_func_search_miss_eig(eig_target, eig_search)
if isempty(eig_target)
    idx_eig_target_loss = [];
    max_mismatch = 0;
else
    flg_found = zeros(length(eig_target),1);
    max_mismatch = 0;
    for ii = 1:length(eig_target)
        mismatch = min(abs(eig_target(ii) - eig_search));
        max_mismatch = max([max_mismatch, mismatch]);
        if any(abs(eig_target(ii) - eig_search) < 1e-4)
            flg_found(ii) = 1;
        end
    end
    idx_eig_target_loss = find(flg_found == 0);
end