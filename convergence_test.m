%% convergence_test.m
% Verify that within_host_model_2nd_order achieves O(h^2) convergence while
% the original first-order scheme achieves O(h^1), using Richardson
% extrapolation with a fine-grid reference solution.
%
% Runs both solvers at multiple step sizes and plots log(error) vs log(h).
% The expected slopes are 1 (first-order) and 2 (second-order).

global P
baseline_parameter_set;

%% Reference solution at very fine step size
h_ref = 0.0625;  % reference step (half of the standard step)

X_max   = 200*24;   % shorter run for speed during convergence test
tau_max = 20*24;
x_ref   = (0:h_ref:X_max)';
nx_ref  = length(x_ref);
tau_ref = (0:h_ref:tau_max)';
ntau_ref = length(tau_ref);

B0_ref = P.Bstar;
M0_ref = 0;
I0_ref = ones(1, ntau_ref);
I0_ref(floor(48/h_ref)+1:end) = 0;
initial_innoc = 0.06;
I0_ref = initial_innoc * I0_ref / (h_ref * trapz(I0_ref));
IG0_ref = zeros(1, ntau_ref);
G0_ref  = 0;
A0_ref  = 0;
CC_ref  = P.c * ones(1, nx_ref);

disp('Computing reference solution (h = 0.0625)...');
[~, ~, G_ref, ~, ~] = within_host_model_2nd_order( ...
    h_ref, 0, X_max, tau_max, B0_ref, M0_ref, I0_ref, IG0_ref, G0_ref, A0_ref, CC_ref);

%% Test step sizes
h_test = [2, 1, 0.5, 0.25, 0.125];

err_1st = zeros(size(h_test));
err_2nd = zeros(size(h_test));

for k = 1:length(h_test)
    h = h_test(k);
    x   = (0:h:X_max)';
    nx  = length(x);
    tau = (0:h:tau_max)';
    ntau = length(tau);

    B0 = P.Bstar;
    M0 = 0;
    I0 = ones(1, ntau);
    I0(floor(48/h)+1:end) = 0;
    I0 = initial_innoc * I0 / (h * trapz(I0));
    IG0 = zeros(1, ntau);
    G0  = 0;
    A0  = 0;
    CC  = P.c * ones(1, nx);

    fprintf('h = %.4f hr: ', h);

    % --- 1st-order scheme (original) ---
    [~, ~, ~, ~, ~, ~] = deal([]);  % clear previous outputs
    [B1, M1, I1, IG1, G1_sol, A1] = within_host_model( ...
        h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);

    % --- 2nd-order scheme (new) ---
    [~, ~, G2_sol, ~, ~] = within_host_model_2nd_order( ...
        h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);

    % Interpolate reference solution to coarse grid for error comparison
    G_ref_interp = interp1(x_ref, G_ref, x, 'linear', 'extrap');

    % L-infinity error in G
    err_1st(k) = max(abs(G1_sol - G_ref_interp));
    err_2nd(k) = max(abs(G2_sol - G_ref_interp));

    fprintf('err_1st = %.4e,  err_2nd = %.4e\n', err_1st(k), err_2nd(k));
end

%% Compute and display convergence rates
fprintf('\nConvergence rates (slope in log-log plot):\n');
rates_1st = diff(log(err_1st)) ./ diff(log(h_test));
rates_2nd = diff(log(err_2nd)) ./ diff(log(h_test));
for k = 1:length(h_test)-1
    fprintf('  h: %.3f -> %.3f:  1st-order rate = %.2f,  2nd-order rate = %.2f\n', ...
        h_test(k), h_test(k+1), rates_1st(k), rates_2nd(k));
end

%% Plot
figure;
loglog(h_test, err_1st, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
loglog(h_test, err_2nd, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);

% Reference lines showing O(h) and O(h^2) slopes
h_line = logspace(log10(min(h_test)), log10(max(h_test)), 50);
scale1 = err_1st(end) / h_test(end)^1;
scale2 = err_2nd(end) / h_test(end)^2;
loglog(h_line, scale1 * h_line.^1, 'b--', 'LineWidth', 1);
loglog(h_line, scale2 * h_line.^2, 'r--', 'LineWidth', 1);

xlabel('Step size h (hours)');
ylabel('Max error in G(x) vs reference');
legend('1st-order scheme', '2nd-order scheme', 'O(h)', 'O(h^2)', 'Location', 'northwest');
title('Convergence test: G(x) error vs step size');
set(gca, 'TickDir', 'out');
box off;
grid on;

disp('Done. Slope ≈ 1 for first-order, ≈ 2 for second-order.');
