%% speed_test.m
% Compare wall-clock time of the original first-order solver vs the
% second-order solver across a range of step sizes.
%
% For a fair comparison, both solvers are run with the same h and the same
% initial conditions. The 2nd-order solver does ~2x more work per step but
% uses only O(ntau) memory instead of O(nx*ntau), so cache effects favour
% it especially at small h.

global P
baseline_parameter_set;

%% Shared setup
X_max   = 1000*24;
tau_max = 20*24;

h_values  = [0.5, 0.25, 0.125, 0.0625];  % step sizes (hours)
n_repeats = 3;                    % repeats per (solver, h) pair to reduce timing noise

t_1st = zeros(length(h_values), n_repeats);
t_2nd = zeros(length(h_values), n_repeats);

fprintf('%-8s  %-12s  %-12s  %-10s\n', 'h (hr)', 't_1st (s)', 't_2nd (s)', 'speedup');
fprintf('%s\n', repmat('-', 1, 48));

for k = 1:length(h_values)
    h = h_values(k);

    tau  = (0:h:tau_max)';
    ntau = length(tau);
    x    = (0:h:X_max)';
    nx   = length(x);

    % Initial conditions
    B0  = P.Bstar;
    M0  = 0;
    I0  = ones(1, ntau);
    I0(floor(48/h)+1:end) = 0;
    I0  = 0.06 * I0 / (h * trapz(I0));
    IG0 = zeros(1, ntau);
    G0  = 0;
    A0  = 0;
    CC  = P.c * ones(1, nx);

    % Time 1st-order solver
    for r = 1:n_repeats
        tic;
        within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);
        t_1st(k, r) = toc;
    end

    % Time 2nd-order solver
    for r = 1:n_repeats
        tic;
        within_host_model_2nd_order(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);
        t_2nd(k, r) = toc;
    end

    med1 = median(t_1st(k,:));
    med2 = median(t_2nd(k,:));
    fprintf('%-8.4f  %-12.3f  %-12.3f  %-10.2fx\n', h, med1, med2, med1/med2);
end

fprintf('\n');
fprintf('Note: at equal h, the 2nd-order solver does ~2x more work per step.\n');
fprintf('Its advantage comes from (a) smaller h needed for same accuracy and\n');
fprintf('(b) rolling-window memory (O(ntau) vs O(nx*ntau)) improving cache use.\n');

%% Plot median times
figure;
semilogy(h_values, median(t_1st, 2), 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(h_values, median(t_2nd, 2), 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Step size h (hours)');
ylabel('Wall-clock time (s)');
legend('1st-order solver', '2nd-order solver', 'Location', 'northwest');
title('Speed comparison: old vs new solver (1000-day simulation)');
set(gca, 'XDir', 'reverse');   % smaller h (more expensive) on right
set(gca, 'TickDir', 'out');
box off;
grid on;
