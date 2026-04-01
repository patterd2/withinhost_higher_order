%% run_strategy_optimization_v2.m
% Redesigned strategy optimiser for nonconstant gametocyte investment.
%
% Improvements over run_strategy_optimization.m:
%   1. Persistent basis-matrix cache  — importdata called once per worker,
%      not once per objective evaluation (biggest single speedup).
%   2. Parallel seed sweep (parfor)   — seeds are independent; wall-clock
%      time scales with 1/num_workers.
%   3. Latin Hypercube initial seeds  — better space-filling than pure
%      random (requires Statistics & ML Toolbox; falls back to rand).
%   4. Two-phase optimisation         — coarse parallel sweep (MaxIter=80)
%      identifies the best seeds cheaply; a refined serial pass
%      (fmincon if available, else fminsearch with tighter tolerances)
%      polishes the top-K survivors to high accuracy.
%
% NB: baseline_parameter_set must NOT be called inside
%     withinhost_model_optimization — it is called once here and
%     replicated to workers via the P_local closure variable.

tic

%% Graphics defaults
set(0,'defaultTextFontName', 'Arial');
set(0,'defaultaxesfontsize', 20);
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%% Parameters (must match withinhost_model_optimization.m)
global P
baseline_parameter_set;
X_max = 1000*24;
h     = 0.125;      % step size (hr). Use 0.25 for ~2x speedup.
x     = (0:h:X_max)';
nx    = length(x);

%% Optimisation configuration
N_SEEDS      = 60;     % number of random seeds for the coarse phase
N_DIM        = 4;      % spline weight dimension (cubic, no interior knot)
WEIGHT_RANGE = 1.0;    % seeds drawn from [-WEIGHT_RANGE, WEIGHT_RANGE]^N_DIM
K_REFINE     = 5;      % number of top seeds passed to the refinement phase
MAXITER_COARSE  = 80;  % fminsearch iterations in the coarse phase
MAXITER_REFINE  = 500; % iterations in the refinement phase

% Bounds for fmincon (wide — prevents wandering off to infinity)
lb = -5 * ones(1, N_DIM);
ub =  5 * ones(1, N_DIM);

%% -----------------------------------------------------------------------
%% Build initial seed matrix  [N_SEEDS × N_DIM]
%% -----------------------------------------------------------------------
if exist('lhsdesign', 'file') == 2
    % Latin Hypercube Sampling (Statistics & ML Toolbox)
    raw_seeds  = lhsdesign(N_SEEDS, N_DIM);         % in [0,1]^N_DIM
    all_seeds  = WEIGHT_RANGE * (2*raw_seeds - 1);  % rescale to [-W, W]
else
    % Fallback to plain random
    all_seeds  = WEIGHT_RANGE * (rand(N_SEEDS, N_DIM) - 0.5);
end

N_TOTAL = size(all_seeds, 1);

%% -----------------------------------------------------------------------
%% Parallel-pool initialisation
%% -----------------------------------------------------------------------
% P is passed explicitly through a closure so that each parfor worker
% receives a local copy — MATLAB does not share global variables across
% workers, and spmd does not support global declarations.
P_local = P;   % capture struct for parfor closure

use_parallel = ~isempty(ver('parallel'));
if use_parallel
    if isempty(gcp('nocreate'))
        parpool('local');
    end
    n_workers = gcp().NumWorkers;
else
    n_workers = 1;
end

%% -----------------------------------------------------------------------
%% PHASE 1 — Coarse parallel sweep
%% -----------------------------------------------------------------------
fprintf('=== Phase 1: Coarse sweep  (%d seeds, MaxIter=%d, %d workers) ===\n', ...
    N_TOTAL, MAXITER_COARSE, n_workers);

coarse_weights = zeros(N_TOTAL, N_DIM);
coarse_fvals   = zeros(N_TOTAL, 1);

opts_coarse = optimset('Display','off','MaxIter', MAXITER_COARSE);

parfor i = 1:N_TOTAL
    % The closure variable P_local is serialised and sent to each worker.
    % withinhost_model_optimization sets its global P from P_local on entry.
    [w, fv] = fminsearch(@(w) withinhost_model_optimization(w, P_local), ...
                         all_seeds(i,:), opts_coarse);
    coarse_weights(i,:) = w;
    coarse_fvals(i)     = fv;
end

% Sort seeds by fitness (ascending fval = descending fitness)
[coarse_fvals_sorted, sort_idx] = sort(coarse_fvals);
top_seeds = coarse_weights(sort_idx(1:K_REFINE), :);

fprintf('  Top %d seeds after coarse sweep:\n', K_REFINE);
for k = 1:K_REFINE
    fprintf('    [%d]  f = %.5f  weights = [%s]\n', k, -coarse_fvals_sorted(k), ...
        num2str(top_seeds(k,:), '%.4f '));
end
fprintf('\n');

%% -----------------------------------------------------------------------
%% PHASE 2 — Refined local search on top-K survivors
%% -----------------------------------------------------------------------
fprintf('=== Phase 2: Refinement  (top %d seeds, MaxIter=%d) ===\n', ...
    K_REFINE, MAXITER_REFINE);

refine_weights = zeros(K_REFINE, N_DIM);
refine_fvals   = zeros(K_REFINE, 1);

use_fmincon = exist('fmincon', 'file') == 2;   % Optimization Toolbox?
if use_fmincon
    opts_refine = optimoptions('fmincon', ...
        'Display',             'off', ...
        'MaxIterations',       MAXITER_REFINE, ...
        'FunctionTolerance',   1e-7, ...
        'StepTolerance',       1e-7, ...
        'FiniteDifferenceType','central');
    fprintf('  Using fmincon (gradient-based, Optimization Toolbox).\n\n');
else
    opts_refine = optimset('Display','off','MaxIter', MAXITER_REFINE, ...
        'TolFun', 1e-7, 'TolX', 1e-7);
    fprintf('  Optimization Toolbox not found — using fminsearch.\n\n');
end

obj = @(w) withinhost_model_optimization(w, P_local);   % closure, consistent with Phase 1

for k = 1:K_REFINE
    if use_fmincon
        [w, fv] = fmincon(obj, top_seeds(k,:), ...
            [], [], [], [], lb, ub, [], opts_refine);
    else
        [w, fv] = fminsearch(obj, top_seeds(k,:), opts_refine);
    end
    refine_weights(k,:) = w;
    refine_fvals(k)     = fv;
    fprintf('  Seed %d/%d  ->  f = %.5f\n', k, K_REFINE, -fv);
end

%% -----------------------------------------------------------------------
%% Identify the overall best strategy
%% -----------------------------------------------------------------------
[best_fval, best_k] = min(refine_fvals);
best_fitness = -best_fval;
best_weights = refine_weights(best_k, :);

fprintf('\n=== Result ===\n');
fprintf('Optimal fitness:  f = %.5f\n', best_fitness);
fprintf('Optimal weights:  [%s]\n\n', num2str(best_weights, '%.6f '));

elapsed = toc;
fprintf('Total elapsed time: %.1f s  (%.1f min)\n\n', elapsed, elapsed/60);

%% -----------------------------------------------------------------------
%% Save results
%% -----------------------------------------------------------------------
save('optimization_results.mat', ...
    'best_weights', 'best_fitness', ...
    'refine_weights', 'refine_fvals', ...
    'coarse_weights', 'coarse_fvals');
fprintf('Results saved to optimization_results.mat\n');

%% -----------------------------------------------------------------------
%% Reconstruct and plot the optimal strategy
%% -----------------------------------------------------------------------
basis_data = importdata('basisMatrixNoKnots_1000_0.125.txt');
B_matrix   = basis_data.data;   % nx × N_DIM

CC_opt = zeros(nx, 1);
for k = 1:N_DIM
    CC_opt = CC_opt + best_weights(k) * B_matrix(:, k);
end
CC_opt = min(1, max(0, CC_opt));

figure(1);
clf;
plot(x/24, 100*CC_opt, 'Color', [0 0.45 0.70], 'LineWidth', 3);
ylim([0 50]);
xlim([0 600]);
xlabel('infection age x (days)');
ylabel('transmission investment (%)');
xtickformat('%.0f');
ytickformat('percentage');
title('Optimal plastic investment strategy');
set(gca, 'TickDir', 'out');
box off;
