function b = withinhost_model_optimization(spline_weights)
% WITHINHOST_MODEL_OPTIMIZATION  Objective function for strategy optimization.
%   Returns the negative cumulative infectiousness (fitness) for a given
%   set of spline weights defining the investment strategy CC(x).
%   Uses the second-order solver within_host_model_2nd_order.
%
%   Called by fminsearch / fmincon in run_strategy_optimization_v2.m.
%
%   Performance notes
%   -----------------
%   * The basis matrix (192 K rows) and the initial-condition vectors are
%     cached via persistent variables and loaded/computed only once per
%     MATLAB session (or once per parallel worker).
%   * baseline_parameter_set must NOT be called inside this function.
%     The caller is responsible for setting global P before the first call.

global P

%% Persistent cache — built once per session / worker
persistent H_CACHED                    % step size at last cache build
persistent X NX AC PSI                 % grids and no-immunity constants
persistent B0 M0 I0 IG0 G0 A0         % initial conditions (invariant to weights)
persistent BASIS   BASIS_DEG           % basis matrix and its polynomial degree

%% Numerical configuration (must match runner script)
X_max   = 1000*24;
tau_max = 20*24;
h       = 0.125;    % step size in hours (use 0.25 for 2x speedup)

G_threshold = 1;

deg = length(spline_weights);

%% Rebuild grid / IC cache if h has changed (or on first call)
if isempty(H_CACHED) || H_CACHED ~= h

    H_CACHED = h;

    tau  = (0:h:tau_max)';
    ntau = length(tau);
    X    = (0:h:X_max)';
    NX   = length(X);
    AC   = floor(280*24/h) + 1;   % 280-day cap index (no-immunity case)
    PSI  = 1/105;                  % recovery rate (no-immunity case)

    % Initial conditions (same for every strategy)
    B0 = P.Bstar;
    M0 = 0;
    I0 = ones(1, ntau);
    I0(floor(48/h)+1:end) = 0;
    initial_innoc = 0.06;
    I0  = initial_innoc * I0 / (h * trapz(I0));
    IG0 = zeros(1, ntau);
    G0  = 0;
    A0  = 0;

    % Invalidate basis cache when h changes (filenames encode h)
    BASIS     = [];
    BASIS_DEG = -1;

end

%% Lazy-load basis matrix for the requested polynomial degree
% Mapping: number of weights → basis file (each file has that many columns)
%   2 weights → degree-1 (linear) spline
%   3 weights → degree-2 (quadratic) spline
%   4 weights → degree-3 (cubic) spline, no interior knot  [default]
%   5 weights → degree-4 (quartic) spline
if isempty(BASIS) || BASIS_DEG ~= deg
    switch deg
        case 2,  fname = 'basisMatrixNoKnots_degree1_1000_0.125.txt';
        case 3,  fname = 'basisMatrixNoKnots_degree2_1000_0.125.txt';
        case 4,  fname = 'basisMatrixNoKnots_1000_0.125.txt';
        case 5,  fname = 'basisMatrixNoKnots_degree4_1000_0.125.txt';
        otherwise
            error('withinhost_model_optimization: unsupported weight vector length %d ' ...
                  '(scalar = constant strategy; 2–5 = spline basis).', deg);
    end
    tmp       = importdata(fname);
    BASIS     = tmp.data;
    BASIS_DEG = deg;
end

%% Build investment strategy CC from spline weights (length = deg)
if deg == 1
    CC = spline_weights(1) * ones(1, NX);
else
    % Linear combination of basis columns, clipped to [0, 1]
    CC = zeros(NX, 1);
    for k = 1:deg
        CC = CC + spline_weights(k) * BASIS(:, k);
    end
    CC = min(1, max(0, CC))';
end

%% Run second-order solver (only G is needed for fitness)
[~, ~, G, ~, ~] = within_host_model_2nd_order( ...
    h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);

%% Compute fitness (cumulative infectiousness)
length_infection = find(G > G_threshold, 1, 'last');
if isempty(length_infection)
    length_infection = NX;
end

if P.sigma == 0
    % No immunity: fixed-duration fitness with exponential discounting
    b = -simps(X(1:AC), betaHV(G(1:AC)) .* exp(-PSI*X(1:AC)/24)) / 24;
else
    % With immunity: infection ends when G drops below threshold
    b = -simps(X(1:length_infection), betaHV(G(1:length_infection))) / 24;
end
end
