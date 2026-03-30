%% run_strategy_optimization.m
% Script to perform strategy optimization for nonconstant parasite
% investment using the second-order within-host solver.
%
% NB: Ensure parameters and h agree with withinhost_model_optimization.m.

tic

set(0,'defaultTextFontName', 'Arial');
set(0,'defaultaxesfontsize', 20);
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%% Parameters (must match withinhost_model_optimization.m)
baseline_parameter_set;
X_max = 1000*24;
h     = 0.125;     % step size in hours (use 0.25 for 2x speedup)
x     = (0:h:X_max)';

N = 100;  % number of random initial seeds
save_strats    = zeros(N, 4);
max_cum_inf    = zeros(N, 1);

for i = 1:N
    % Random initial seed (4 weights for cubic spline, no knot)
    v = 0.5*(rand(1,4) - 0.5);
    options = optimset('Display', 'off', 'MaxIter', 200);
    [a, funmax] = fminsearch(@withinhost_model_optimization, v, options);
    max_cum_inf(i)         = -funmax;
    save_strats(i, 1:length(a)) = a;
    fprintf('Seed %3d/%3d:  f = %.5f\n', i, N, -funmax);
end
toc

%% Identify and plot the best strategy
[~, opt_strat] = max(max_cum_inf);
a = save_strats(opt_strat, :);

temp1 = importdata('basisMatrixNoKnots_1000_0.125.txt');
CC = min(1, max(0, a(1)*temp1.data(:,1) + a(2)*temp1.data(:,2) + ...
                   a(3)*temp1.data(:,3) + a(4)*temp1.data(:,4)))';

figure(1);
hold on;
plot(x/24, 100*CC, 'LineWidth', 3);
ylim([0 50]);
xlim([0 600]);
xlabel('infection age x (days)');
ylabel('transmission investment (%)');
ytickformat('percentage');
title('Optimal plastic investment strategy');
set(gca, 'TickDir', 'out');
box off;

fprintf('\nOptimal fitness: %.5f\n', max(max_cum_inf));
fprintf('Optimal weights: [%.6f, %.6f, %.6f, %.6f]\n', a(1), a(2), a(3), a(4));
