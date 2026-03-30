%% compare_schemes.m
% Run the original first-order and the new second-order solver with the
% same parameters and plot key state variables side-by-side.
%
% Edit h_old / h_new below to test different step sizes.
% Setting h_new = 0.25 with h_old = 0.125 demonstrates comparable accuracy
% at roughly half the cost.

global P
baseline_parameter_set;

%% Configuration
h_old = 0.125;   % step size for 1st-order solver
h_new = 0.125;   % step size for 2nd-order solver (try 0.25)

X_max   = 1000*24;
tau_max = 20*24;
G_threshold = 1;

%% Initial conditions (old solver)
tau_o  = (0:h_old:tau_max)';
ntau_o = length(tau_o);
nx_o   = length(0:h_old:X_max);

B0_o  = P.Bstar;  M0_o = 0;
I0_o  = ones(1, ntau_o);
I0_o(floor(48/h_old)+1:end) = 0;
I0_o  = 0.06 * I0_o / (h_old * trapz(I0_o));
IG0_o = zeros(1, ntau_o);
CC_o  = P.c * ones(1, nx_o);

%% Initial conditions (new solver)
tau_n  = (0:h_new:tau_max)';
ntau_n = length(tau_n);
nx_n   = length(0:h_new:X_max);

B0_n  = P.Bstar;  M0_n = 0;
I0_n  = ones(1, ntau_n);
I0_n(floor(48/h_new)+1:end) = 0;
I0_n  = 0.06 * I0_n / (h_new * trapz(I0_n));
IG0_n = zeros(1, ntau_n);
CC_n  = P.c * ones(1, nx_n);

%% Run old (1st-order) solver
fprintf('Running 1st-order solver (h = %.4f hr)...\n', h_old);
tic;
[B_old, M_old, I_old, ~, G_old, A_old] = within_host_model( ...
    h_old, 0, X_max, tau_max, B0_o, M0_o, I0_o, IG0_o, 0, 0, CC_o);
t_old = toc;
inf_old = h_old * sum(I_old, 2);   % total iRBC density at each x
fprintf('  Done: %.2f s\n', t_old);

%% Run new (2nd-order) solver
fprintf('Running 2nd-order solver (h = %.4f hr)...\n', h_new);
tic;
[B_new, M_new, G_new, A_new, inf_new] = within_host_model_2nd_order( ...
    h_new, 0, X_max, tau_max, B0_n, M0_n, I0_n, IG0_n, 0, 0, CC_n);
t_new = toc;
fprintf('  Done: %.2f s\n', t_new);

%% Time axes (in days)
x_old = (0:h_old:X_max)' / 24;
x_new = (0:h_new:X_max)' / 24;

%% Infection end points
end_old = find(G_old > G_threshold, 1, 'last');
end_new = find(G_new > G_threshold, 1, 'last');
if isempty(end_old), end_old = length(x_old); end
if isempty(end_new), end_new = length(x_new); end
x_end = max(x_old(end_old), x_new(end_new));

%% Summary diagnostics
fprintf('\nInfection length:   1st-order = %.1f days,   2nd-order = %.1f days\n', ...
    x_old(end_old), x_new(end_new));
fprintf('Cumul. infect.:     1st-order = %.5f,   2nd-order = %.5f\n', ...
    simps(0:h_old:X_max, betaHV(G_old))/24, ...
    simps(0:h_new:X_max, betaHV(G_new))/24);
fprintf('Speedup:            %.2fx\n', t_old / t_new);

%% Interpolate new solver onto old grid for pointwise difference panel
if h_old ~= h_new
    G_new_i   = interp1(x_new, G_new,   x_old, 'linear', 'extrap');
    M_new_i   = interp1(x_new, M_new,   x_old, 'linear', 'extrap');
    inf_new_i = interp1(x_new, inf_new, x_old, 'linear', 'extrap');
else
    G_new_i   = G_new;
    M_new_i   = M_new;
    inf_new_i = inf_new;
end

%% Plotting
set(0, 'defaultaxesfontsize', 13);
figure('Position', [80 80 1400 880]);

lw = 2;
c1 = [0.00, 0.45, 0.70];   % blue  — 1st-order
c2 = [0.85, 0.33, 0.10];   % red   — 2nd-order

% Panel 1: iRBC abundance
subplot(2, 3, 1);
semilogy(x_old, max(inf_old, 1e-10), '-',  'Color', c1, 'LineWidth', lw); hold on;
semilogy(x_new, max(inf_new, 1e-10), '--', 'Color', c2, 'LineWidth', lw);
xlim([0 x_end]);
xlabel('infection age (days)');
ylabel('iRBC abundance');
title('\int I d\tau  (total iRBC)');
legend(sprintf('1st-order  h=%.3g', h_old), sprintf('2nd-order  h=%.3g', h_new), ...
       'Location', 'best');
set(gca, 'TickDir', 'out'); box off;

% Panel 2: Gametocytes
subplot(2, 3, 2);
semilogy(x_old, max(G_old, 1e-10), '-',  'Color', c1, 'LineWidth', lw); hold on;
semilogy(x_new, max(G_new, 1e-10), '--', 'Color', c2, 'LineWidth', lw);
xlim([0 x_end]);
xlabel('infection age (days)');
ylabel('abundance');
title('Gametocytes G(x)');
set(gca, 'TickDir', 'out'); box off;

% Panel 3: Merozoites
subplot(2, 3, 3);
semilogy(x_old, max(M_old, 1e-10), '-',  'Color', c1, 'LineWidth', lw); hold on;
semilogy(x_new, max(M_new, 1e-10), '--', 'Color', c2, 'LineWidth', lw);
xlim([0 x_end]);
xlabel('infection age (days)');
ylabel('abundance');
title('Merozoites M(x)');
set(gca, 'TickDir', 'out'); box off;

% Panel 4: Uninfected RBCs
subplot(2, 3, 4);
plot(x_old, B_old/1e6, '-',  'Color', c1, 'LineWidth', lw); hold on;
plot(x_new, B_new/1e6, '--', 'Color', c2, 'LineWidth', lw);
xlim([0 x_end]);
xlabel('infection age (days)');
ylabel('B  (\times10^6 cells/\muL)');
title('Uninfected RBCs B(x)');
set(gca, 'TickDir', 'out'); box off;

% Panel 5: Immune effector
subplot(2, 3, 5);
plot(x_old, A_old, '-',  'Color', c1, 'LineWidth', lw); hold on;
plot(x_new, A_new, '--', 'Color', c2, 'LineWidth', lw);
xlim([0 x_end]);
xlabel('infection age (days)');
ylabel('A');
title('Immune effector A(x)');
set(gca, 'TickDir', 'out'); box off;

% Panel 6: Pointwise relative difference in G
subplot(2, 3, 6);
rel_diff = abs(G_new_i - G_old) ./ (abs(G_old) + 1e-12);
semilogy(x_old, rel_diff, 'k-', 'LineWidth', lw);
xlim([0 x_end]);
xlabel('infection age (days)');
ylabel('|new - old| / |old|');
title(sprintf('Relative difference in G\n(h_{old}=%.3g, h_{new}=%.3g)', h_old, h_new));
set(gca, 'TickDir', 'out'); box off; grid on;

sgtitle(sprintf('Scheme comparison  |  1st-order: %.2f s  |  2nd-order: %.2f s  |  speedup: %.2fx', ...
    t_old, t_new, t_old/t_new), 'FontSize', 14);
