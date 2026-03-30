function b = withinhost_model_optimization(spline_weights)
% WITHINHOST_MODEL_OPTIMIZATION  Objective function for strategy optimization.
%   Returns the negative cumulative infectiousness (fitness) for a given
%   set of spline weights defining the investment strategy CC(x).
%   Uses the second-order solver within_host_model_2nd_order.
%
%   Called by fminsearch in run_strategy_optimization.m.

global P

%% Numerical configuration (must match run_strategy_optimization.m)
X_max   = 1000*24;
tau_max = 20*24;
h       = 0.125;     % step size — use 0.25 for 2x speedup with comparable accuracy

tau  = (0:h:tau_max)';
ntau = length(tau);
x    = (0:h:X_max)';
nx   = length(x);
ac   = floor(280*24/h) + 1;   % index for 280-day cap (no-immunity case)
psi  = 1/105;                  % constant recovery rate (no-immunity case)

G_threshold = 1;

baseline_parameter_set;

%% Initial conditions
B0 = P.Bstar;
M0 = 0;
I0 = ones(1, ntau);
I0(floor(48/h)+1:end) = 0;
initial_innoc = 0.06;
I0 = initial_innoc * I0 / (h * trapz(I0));
IG0 = zeros(1, ntau);
G0  = 0;
A0  = 0;

%% Build investment strategy from spline weights
if isscalar(spline_weights)
    CC = spline_weights * ones(1, nx);
elseif length(spline_weights) == 2
    temp1 = importdata('basisMatrixNoKnots_degree1_1000_0.125.txt');
    CC = min(1, max(0, spline_weights(1)*temp1.data(:,1) + ...
                       spline_weights(2)*temp1.data(:,2)));
elseif length(spline_weights) == 3
    temp1 = importdata('basisMatrixNoKnots_degree2_1000_0.125.txt');
    CC = min(1, max(0, spline_weights(1)*temp1.data(:,1) + ...
                       spline_weights(2)*temp1.data(:,2) + ...
                       spline_weights(3)*temp1.data(:,3)));
elseif length(spline_weights) == 4
    temp1 = importdata('basisMatrixNoKnots_1000_0.125.txt');
    CC = min(1, max(0, spline_weights(1)*temp1.data(:,1) + ...
                       spline_weights(2)*temp1.data(:,2) + ...
                       spline_weights(3)*temp1.data(:,3) + ...
                       spline_weights(4)*temp1.data(:,4)));
elseif length(spline_weights) == 5
    temp1 = importdata('basisMatrixNoKnots_degree4_1000_0.125.txt');
    CC = min(1, max(0, spline_weights(1)*temp1.data(:,1) + ...
                       spline_weights(2)*temp1.data(:,2) + ...
                       spline_weights(3)*temp1.data(:,3) + ...
                       spline_weights(4)*temp1.data(:,4) + ...
                       spline_weights(5)*temp1.data(:,5)));
end

%% Run second-order solver (only G is needed for fitness)
[~, ~, G, ~, ~] = within_host_model_2nd_order( ...
    h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);

%% Compute fitness (cumulative infectiousness)
length_infection = find(G > G_threshold, 1, 'last');
if isempty(length_infection)
    length_infection = length(x);
end

if P.sigma == 0
    % No immunity: fixed-duration fitness
    b = -simps(x(1:ac), betaHV(G(1:ac)) .* exp(-psi*x(1:ac)/24)) / 24;
else
    % With immunity: infection ends when G drops below threshold
    b = -simps(x(1:length_infection), betaHV(G(1:length_infection))) / 24;
end
end
