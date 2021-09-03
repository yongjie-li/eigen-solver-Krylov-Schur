load mat_ABCD_case39.mat

% shift point for Shift-invert transform
x_shift = 0;
y_shift = 1.5;
r_search = 2; % the radius of circle with center at (x_shitf, y_shift), within which we search eigenvalues

% setting of algorithm pars.
opt_solver.dim_subspace = 80;
opt_solver.n_restart_max = 15;
opt_solver.tol_residue = 1e-6;

% search eigenvalue using Krylov-Schur 
idx_eig_target = abs(e_base - (x_shift + 1j * y_shift)) < r_search; % comparing eigenvalues within search circle
mat_ABCD.eig_target = e_base(idx_eig_target); % field 'eig_target' is not necessary

[v, e_converge, e_approximate, r_cover, flg_eig_dist] = func_search_eig_Krylov_Schur(x_shift, y_shift, r_search, mat_ABCD, opt_solver);

% plot result
id_figure = figure();
plot(real(e_base), imag(e_base), 's');
hold on
plot(real(e_converge), imag(e_converge), 'x');
plot(x_shift, y_shift, 'o');
aux_func_plot_circle(x_shift, y_shift, r_cover, id_figure, '--')
legend('Accurate eigenvalue', 'Found eigenvalue by KS', 'Shift point', 'Searched circle');
xlim([x_shift - r_search - 1, x_shift + r_search + 1]);
ylim([y_shift - r_search - 1, y_shift + r_search + 1]);
hold off
