function [h_mapped, v_subspace, w, r, n_converge] = func_Arnoldi_factorize(h_mapped, v_subspace, dim_subspace, idx_v_start, lu_shift_A, lu_eqn1_coef, mat_state_B, mat_state_C, tol_converge)
    n_converge = 0;
    for idx_v = idx_v_start : dim_subspace
        w = func_multiply_shifted_ABCD(lu_shift_A, lu_eqn1_coef, mat_state_B, mat_state_C, v_subspace(:, idx_v));
        h_mapped(1:idx_v, idx_v) = v_subspace(:, 1:idx_v)' * w;
        w = w - v_subspace(:, 1:idx_v) * h_mapped(1:idx_v, idx_v);
        r = norm(w, 2);
        w = w / r;
        if r > tol_converge % continue Arnoldi factorize
            if idx_v < dim_subspace
                h_mapped(idx_v + 1,idx_v) = r;
                v_subspace(:, idx_v + 1) = w;               
            end
        else % converge
            n_converge = idx_v;
            disp('Arnoldi procedure terminates due to zero residue');
            break;
        end
    end
end