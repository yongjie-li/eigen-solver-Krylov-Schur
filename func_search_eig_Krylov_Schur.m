function [v, e_converge, e_approximate, r_cover, flg_eig_dist] = func_search_eig_Krylov_Schur(x_center_circle, y_center_circle, r_search_circle, mat_ABCD, opt_solver)
    %% input pars:
    % x_center_circle: real part of shift point
    % y_center_circle: imag. part of shift point
    % r_search_circle: radius of the circle to be searched
    
    % mat_ABCD is struct with fields:
    % .A 
    % .B
    % .C
    % .D
    % [dx/dt] = [ A B ] [x]
    % [ 0   ] = [ C D ] [y]
    % attention: A, B, C, D are with sparse matrix format
    % opt_solver is struct with fields:
    % .dim_subspace : dimension of Krylov subspace
    % .n_restart_max : max number of restarts
    % .tol_residue : tolerance for Ritz value convergence
    %% output
    % v: eigen-vectors of converged eigenvalues
    % e_converge: converged Ritz value
    % e_approximate: all Ritz values, including un-converged ones
    % r_cover: radius of circles, within which all eigenvalues are found
    % flg_eig_dist: return status of searching eigenvalue: 1: search circle is covered 2: kept subspace exceeds whole subspace; 3: too many void restart; 4: shift point too close to eig 
    
    flg_eig_dist = 0;  % 1: search circle is covered 2: kept subspace exceeds whole subspace; 3: too many void restart; 4: shift point too close to eig 
    
    %% build data for shifted matrix multiplication
    % eqn1. [C*(A-sI)^-1*B - D] * w = C*(A-sI)^-1*v
    % eqn2. u = (A-sI)^-1*(v - B*w)
    
    point_shift = x_center_circle + 1j * y_center_circle;
    dim_A = size(mat_ABCD.A,1);
    A = mat_ABCD.A;
    B = mat_ABCD.B;
    C = mat_ABCD.C;
    D = mat_ABCD.D;
    
    if isfield(mat_ABCD, 'eig_target')
        eig_target = mat_ABCD.eig_target;
        flg_have_accurate_eig = 1;
    else
        flg_have_accurate_eig = 0;
    end
    
    lu_shift_A = klu(A - speye(dim_A) * point_shift);
    lu_eqn1_coef = klu(sparse(C * klu(lu_shift_A, '\', full(B)) - D));
    %% krylov-schur iteration
    % init
    dim_subspace = opt_solver.dim_subspace;
    n_restart_max = opt_solver.n_restart_max;
    tol_residue = opt_solver.tol_residue;
    v_subspace = zeros(dim_A, dim_subspace);
    h_mapped = zeros(dim_subspace, dim_subspace);
    
    v_subspace(:,1) = rand(dim_A,1);
    v_subspace(:,1) = v_subspace(:,1)/norm(v_subspace(:,1),2);
    
    % Krylov-Schur iteration starts
    cnt_restart = 0;
    n_Ritz_guard = ceil(dim_subspace/10);
    n_Ritz_max_target = ceil(dim_subspace / 10 * 8);
    dim_subspace_kept = n_Ritz_max_target + n_Ritz_guard;
    n_converged_Ritz = 0;
    r_cover = 0;
    cnt_restart_void = 0;
    n_restart_void_max_without_converge = 2;
    n_restart_void_max_with_converge = 4;
    tol_distance_eig_shift_point = 1e-4 - 1e-10;
    scalor_safe_r_cover = 1 - 0.05 * n_restart_void_max_with_converge;
    while(1)
        if cnt_restart == 0 % first Arnoldi factorize
            idx_v_start = 1; % otherwise, idx_v_start is determined by the dimension of kept subspace
        end
        [h_mapped, v_subspace, w, r, n_converged_Ritz_Arnoldi] = func_Arnoldi_factorize(h_mapped, v_subspace, dim_subspace, idx_v_start, lu_shift_A, lu_eqn1_coef, B, C, tol_residue);
        if r < tol_residue % Krylov-Schur iteration terminate due to subspace is consumed
            [vec_h_mapped, eig_h_mapped] = eig(h_mapped(1:n_converged_Ritz_Arnoldi, 1:n_converged_Ritz_Arnoldi));
            e_converge = 1./diag(eig_h_mapped) + point_shift;
            e_approximate = [];
            v = v_subspace(:,1:n_converged_Ritz_Arnoldi)*vec_h_mapped;
            r_cover = max(abs(e_converge - point_shift));
            break;
        else % process the Arnoldi factorization
            % find converged Ritz values
            [Q_schur, h_mapped] = func_part_Schur_factorize(h_mapped, n_converged_Ritz, dim_subspace - n_converged_Ritz);
            v_subspace = v_subspace * Q_schur;
            r_b = zeros(1, dim_subspace); % residual vector
            r_b(1, dim_subspace) = r;
            r_b = r_b * Q_schur;
            
            n_converged_Ritz_previous = n_converged_Ritz;
            n_converged_Ritz  = 0;
            for ii = 1:dim_subspace % find converged Ritz values
                if abs(sum(r_b(ii))) < tol_residue
                    n_converged_Ritz = n_converged_Ritz + 1;
                else
                    break;
                end
            end
            
            [r_covered_new, id_cluster_Ritz, n_Ritz_max_target] = func_process_Ritz(diag(h_mapped), n_converged_Ritz, n_Ritz_max_target, dim_subspace_kept);
%             r_cover = max([r_covered_new, r_cover]);
            r_cover = r_covered_new;
            dim_subspace_kept_old = dim_subspace_kept;
            dim_subspace_kept = n_Ritz_max_target + n_Ritz_guard;
            %
            if n_converged_Ritz == n_converged_Ritz_previous && dim_subspace_kept == dim_subspace_kept_old
                cnt_restart_void = cnt_restart_void + 1;
            else
                cnt_restart_void = 0; % reset
            end
            
            Ritz_converged = diag(h_mapped(1:n_converged_Ritz,1:n_converged_Ritz));
            Ritz_valid = diag(h_mapped(1:n_Ritz_max_target,1:n_Ritz_max_target));
            if any(1./abs(Ritz_converged) <= tol_distance_eig_shift_point)
                flg_big_Ritz = 1./abs(Ritz_converged) < tol_distance_eig_shift_point;
                r_cover = max(1./abs(Ritz_converged(flg_big_Ritz))) + 1e-4;
                e_approximate = []; % approximated Ritz values are not reliable
                flg_terminate = 1;
                flg_eig_dist = 4;
            elseif any(1./abs(Ritz_valid) <= tol_distance_eig_shift_point)
                flg_big_Ritz = 1./abs(Ritz_valid) < tol_distance_eig_shift_point;
                r_cover = max(1./abs(Ritz_valid(flg_big_Ritz))) + 1e-4;
                e_approximate = []; % approximated Ritz values are not reliable
                flg_terminate = 1;
                flg_eig_dist = 4;
            elseif r_cover >= r_search_circle || dim_subspace_kept >= dim_subspace || cnt_restart >= n_restart_max || (cnt_restart_void >= n_restart_void_max_without_converge && n_converged_Ritz == 0) || cnt_restart_void >=  n_restart_void_max_with_converge % Krylov-Schur iteration terminates
                Ritz_approximate = diag(h_mapped((n_converged_Ritz + 1):n_Ritz_max_target, (n_converged_Ritz + 1):n_Ritz_max_target));
                e_approximate = 1./Ritz_approximate + point_shift;
                if dim_subspace_kept >= dim_subspace
                    disp('!kept subspace expands to whole subspace');
                end
                flg_terminate = 1;
                if r_cover >= r_search_circle
                    flg_eig_dist = 1;
                elseif dim_subspace_kept >= dim_subspace
                    flg_eig_dist = 2;
                elseif (cnt_restart_void >= n_restart_void_max_without_converge && n_converged_Ritz == 0) || cnt_restart_void >=  n_restart_void_max_with_converge
                    flg_eig_dist = 3;
                end
            else
                flg_terminate = 0;
            end
            if flg_terminate == 1
                [v_schur, eig_schur] = eig(h_mapped(1:n_converged_Ritz,1:n_converged_Ritz));
                e_converge = 1./diag(eig_schur) + point_shift;
                v = v_subspace(:, 1:n_converged_Ritz) * v_schur;
                r_cover = r_cover * scalor_safe_r_cover;
                if flg_have_accurate_eig == 1 && r_cover > 1e-5
                    flg_in_search_circle = func_is_covered_by_circles(real(eig_target), imag(eig_target), [x_center_circle, y_center_circle, r_cover]);
                    eig_cmp = eig_target(flg_in_search_circle);
                    [idx_miss_eig, max_mismatch] = aux_func_search_miss_eig(eig_cmp, e_converge);
                    if ~isempty(idx_miss_eig)
                        disp(['found missed eigenvalue with max mismatch: ', num2str(max_mismatch)]);
                        figure(3);
                        hold off
                        plot(real(eig_target), imag(eig_target), 'o');
                        hold on
                        plot(real(eig_target(idx_miss_eig)), imag(eig_target(idx_miss_eig)), 's');
                        plot(real(e_converge), imag(e_converge), 'x');
                        plot(real(e_approximate), imag(e_approximate), 'd');
                        aux_func_plot_circle(x_center_circle, y_center_circle, r_cover, 3, 'k--');
                        hold off
                    end
                end
                break;
            else
                [Q_reorder, h_mapped] = ordschur(eye(dim_subspace), h_mapped, id_cluster_Ritz);
                v_subspace = v_subspace * Q_reorder;
                r_b = r_b * Q_reorder;

                idx_v_start = dim_subspace_kept + 1;
                h_mapped(idx_v_start,:) = r_b;
                v_subspace(:, idx_v_start) = w;
                cnt_restart = cnt_restart + 1;
            end
        end
    end
    % evaluate coverage of Ritz values
    if r_cover == 0
        if ~isempty(e_approximate)
            r_cover = min(abs(e_approximate - point_shift)) * scalor_safe_r_cover;
        end
    end
%     if r_cover < tol_residue
%         r_cover = tol_residue;
%     end
    assert(r_cover > 0);
    disp([num2str(length(e_converge)),' eigenvalues converge with reliable searched radius: ',num2str(r_cover),' in ',num2str(cnt_restart),' restarts']);
end