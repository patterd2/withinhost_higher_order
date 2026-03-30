%% run_full_model.m
% Main script to run a single simulation of the within-host model using the
% second-order predictor-corrector solver.

tic
global P
set(0,'defaultTextFontName', 'Arial');
set(0,'defaultaxesfontsize', 20);
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%% Choose constant or nonconstant investment strategy
RUN_constant    = 1;
RUN_nonconstant = 0;
RUN_degree      = 3;  % degree of the polynomial spline (1,2,3,4 - 3 baseline)

%% Numerical configuration
X_max   = 1000*24;   % max infection age (hours)
tau_max = 20*24;     % max age-of-infection tracked (hours)
G_threshold = 1;     % gametocyte threshold to define end of infection

% Step size: can use h = 0.25 with 2nd-order scheme for same accuracy as
% h = 0.125 with the original 1st-order scheme.
h = 0.125;           % time/age step size (hours) — change to 0.25 for 2x speedup

x    = (0:h:X_max)';
nx   = length(x);
tau  = (0:h:tau_max)';
ntau = length(tau);

% Set model parameters
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

%% Set the parasite investment strategy
if RUN_constant
    CC = P.c * ones(1, nx);
else
    if RUN_degree == 3
        temp1 = importdata('basisMatrixNoKnots_1000_0.125.txt');
        CC1 = temp1.data(:,1);
        CC2 = temp1.data(:,2);
        CC3 = temp1.data(:,3);
        CC4 = temp1.data(:,4);
        w1 = 0.167061351199006;
        w2 = 0.401066392028436;
        w3 = -1.555469632674011;
        w4 = 4.000060421096212;
        CC = min(1, max(0, w1*CC1 + w2*CC2 + w3*CC3 + w4*CC4))';
    else
        CC = P.c * ones(1, nx);
    end
end

%% Run the second-order solver
[B, M, G, A, infection_level] = within_host_model_2nd_order( ...
    h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);

%% Plotting and diagnostics
standard_plotting_2nd_order;

disp(['Cumulative infectiousness (no threshold): ', ...
      num2str(simps(0:h:X_max, betaHV(G))/24)]);
toc
