%% 1D_fitness_landscape.m
% Compute the fitness landscape f(c) for constant transmission investment
% strategies c in [0, 0.6] using the second-order within-host solver.
%
% Fitness f is cumulative host infectiousness integrated over the lifespan
% of infection:
%   - With immunity (P.sigma > 0): integral of betaHV(G) from x=0 to the
%     last time G > G_threshold
%   - Without immunity (P.sigma == 0): integral of betaHV(G)*exp(-psi*x/24)
%     up to a 280-day cap, with psi = 1/105 (constant recovery rate)
%
% Uses h = 0.125 hr (change to 0.25 for ~2x speedup with similar accuracy).
% No large I/IG matrices are stored; each solver call uses O(ntau) memory.

tic
global P
baseline_parameter_set;

set(0, 'defaultTextFontName',           'Arial');
set(0, 'defaultaxesfontsize',           18);
set(0, 'defaultAxesTickLabelInterpreter','none');
set(0, 'defaulttextinterpreter',        'none');
set(0, 'defaultAxesXGrid',              'off');
set(0, 'defaultAxesYGrid',              'off');
set(0, 'defaultAxesTickDir',            'out');
set(0, 'defaultAxesLineWidth',          1.5);

%% Configuration
h           = 0.125;       % step size (hr) — use 0.25 for ~2x speedup
X_max       = 1000*24;     % max infection age (hr)
tau_max     = 20*24;       % max age-of-infection tracked (hr)
G_threshold = 1;           % gametocyte threshold defining end of infection

% c sweep: 61 values from 0.01 to 0.60
% (use 0:0.005:0.6 for finer resolution at ~2x cost)
c_values = 0.01:0.01:0.6;
nc       = length(c_values);

% No-immunity parameters (used when P.sigma == 0)
psi = 1/105;               % constant host recovery rate (days^-1)
ac  = floor(280*24/h) + 1; % index cap at day 280

%% Shared grids and initial conditions (fixed across all runs)
x    = (0:h:X_max)';
nx   = length(x);
tau  = (0:h:tau_max)';
ntau = length(tau);

B0  = P.Bstar;
M0  = 0;
I0  = ones(1, ntau);
I0(floor(48/h)+1:end) = 0;
I0  = 0.06 * I0 / (h * trapz(I0));
IG0 = zeros(1, ntau);
G0  = 0;
A0  = 0;

%% Pre-allocate result arrays
cum_inf  = zeros(1, nc);   % fitness f(c)
duration = zeros(1, nc);   % infection duration (days)

%% Sweep over c values
fprintf('Computing fitness landscape: %d values of c in [%.2f, %.2f]\n', ...
    nc, c_values(1), c_values(end));
fprintf('h = %.4f hr,  X_max = %d days\n\n', h, X_max/24);

for ii = 1:nc

    P.c = c_values(ii);
    CC  = P.c * ones(1, nx);

    [~, ~, G, ~, ~] = within_host_model_2nd_order( ...
        h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);

    % Find infection end point
    rec_time = find(G > G_threshold, 1, 'last');
    if isempty(rec_time)
        rec_time = 1;
    end
    duration(ii) = x(rec_time) / 24;   % convert hr -> days

    % Compute fitness
    if P.sigma > 0
        % With adaptive immunity: integrate over actual infection lifespan
        cum_inf(ii) = simps(x(1:rec_time), betaHV(G(1:rec_time))) / 24;
    else
        % No immunity: fixed-duration integral with exponential recovery discount
        cum_inf(ii) = simps(x(1:ac), betaHV(G(1:ac)) .* exp(-psi*x(1:ac)/24)) / 24;
    end

    if mod(ii, round(nc/10)) == 0
        fprintf('  %3.0f%%  c = %.3f  f = %.4f  duration = %.1f days\n', ...
            100*ii/nc, P.c, cum_inf(ii), duration(ii));
    end
end

% Restore baseline c
baseline_parameter_set;

fprintf('\nDone. Total time: %.1f s\n', toc);

%% Identify optimal strategy
[f_max, idx_opt] = max(cum_inf);
c_opt = c_values(idx_opt);
fprintf('Optimal constant strategy:  c* = %.4f (%.2f%%)  ->  f = %.5f\n', ...
    c_opt, 100*c_opt, f_max);

%% Plotting
invest_pct = 100 * c_values;   % percentage for x-axis

figure('Position', [100 100 1000 450]);

% Panel 1: Fitness vs c
subplot(1, 2, 1);
plot(invest_pct, cum_inf, '-', 'Color', [0.00 0.45 0.70], 'LineWidth', 3);
hold on;
scatter(100*c_opt, f_max, 120, 'k', 'filled');
xlim([0 max(invest_pct)]);
ylim([0 1.1*f_max]);
xtickformat('percentage');
xlabel('constant transmission investment c');
ylabel('cumulative infectiousness f');
LimitsX = xlim; LimitsY = ylim;
title('Fitness landscape', 'FontWeight', 'normal', ...
      'HorizontalAlignment', 'left', 'Position', [LimitsX(1), LimitsY(2)]);
set(gca, 'TickDir', 'out');
box off;

% Panel 2: Infection duration vs c
subplot(1, 2, 2);
plot(invest_pct, duration, '--', 'Color', [0.47 0.67 0.19], 'LineWidth', 3);
hold on;
scatter(100*c_opt, duration(idx_opt), 120, 'k', 'filled');
xlim([0 max(invest_pct)]);
xtickformat('percentage');
xlabel('constant transmission investment c');
ylabel('infection duration (days)');
LimitsX = xlim; LimitsY = ylim;
title('Infection duration', 'FontWeight', 'normal', ...
      'HorizontalAlignment', 'left', 'Position', [LimitsX(1), LimitsY(2)]);
set(gca, 'TickDir', 'out');
box off;
