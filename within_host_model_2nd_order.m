function [B, M, G, A, infection_level] = within_host_model_2nd_order(h, t0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC)
% WITHIN_HOST_MODEL_2ND_ORDER  Second-order predictor-corrector solver for the
%   age-structured within-host malaria model.
%
%   Uses a Crank-Nicolson (trapezoidal) scheme along characteristics for
%   the age-structured PDE variables I and IG, combined with a semi-implicit
%   Crank-Nicolson update for the ODE variables B, M, G, A.  The nonlinear
%   cross-coupling terms are handled via a first-order predictor step
%   (identical to the original scheme) followed by a corrector step that
%   averages source terms at step n and the predicted step n+1.
%
%   Memory: only two rows (current and predicted) of I and IG are stored at
%   any time, replacing the O(nx*ntau) storage of the original solver.
%
%   Inputs:
%     h         - time/age step size (hrs), same for both time scales
%     t0        - initial time (usually 0)
%     X_max     - maximum infection age (hrs)
%     tau_max   - maximum age-of-infection tracked (hrs)
%     B0,M0     - initial scalar values for B, M
%     I0        - initial row vector I(t0, tau), length ntau
%     IG0       - initial row vector IG(t0, tau), length ntau
%     G0,A0     - initial scalar values for G, A
%     CC        - investment strategy vector, length nx
%
%   Outputs:
%     B, M, G, A        - solution vectors over x (length nx)
%     infection_level   - h * sum_tau I(x,tau), total iRBC density (length nx)

global P

tau = (0:h:tau_max)';
ntau = length(tau);
x = (t0:h:X_max)';
nx = length(x);

%% Allocate output arrays
% Initialise with zeros: entries beyond the early-exit point remain zero,
% which is biologically correct (infection has ended).
B = zeros(nx,1);
M = zeros(nx,1);
G = zeros(nx,1);
A = zeros(nx,1);
infection_level = zeros(nx,1);

%% Assign initial conditions
B(1) = B0;
M(1) = M0;
G(1) = G0;
A(1) = A0;

I_curr  = I0(:)';   % current row of I,  length ntau
IG_curr = IG0(:)';  % current row of IG, length ntau
infection_level(1) = h * sum(I_curr);

%% Precalculate age-dependent hazard rate vectors (computed once, used every step)
% Gamma1(j): burst hazard at age tau(j), used for current-row integrals
Gamma1 = gamma_fun(tau, h)';             % 1 x ntau

% Gamma2(j): burst hazard at age tau(j+1), i.e. shifted by one, for I corrector denominator
Gamma2 = gamma_fun(tau(2:end)', h);      % 1 x (ntau-1)

% GammaG1(j): gametocyte maturation hazard at age tau(j)
GammaG1 = gamma_G(tau, h)';             % 1 x ntau

% GammaG2(j): gametocyte maturation hazard at age tau(j+1)
GammaG2 = gamma_G(tau(2:end)', h);      % 1 x (ntau-1)

%% Early-exit parameters
% Once G has peaked and dropped back below G_stop, the infection is
% effectively over and all further fitness contributions are negligible.
% Remaining output entries stay zero (set by initialisation above).
G_stop   = 0.01;            % exit threshold for G
n_warmup = floor(10*24/h);  % don't exit before day 10 (G still rising)
G_peaked = false;           % becomes true once G first exceeds G_stop

%% Main time-stepping loop
for n = 1:nx-1

    % ------------------------------------------------------------------
    % Precompute shared integrals at step n (used in both predictor and
    % corrector stages)
    % ------------------------------------------------------------------
    bI_n   = h * sum(Gamma1  .* I_curr);   % beta integral: sum_j Gamma1(j)*I(n,j)*h
    bIG_n  = h * sum(GammaG1 .* IG_curr);  % gametocyte integral
    sumI_n = h * sum(I_curr);              % total iRBC load (for A equation)

    sigma_A_n = P.sigma * (1 - exp(-P.theta * A(n)));  % immune removal rate at step n

    % ------------------------------------------------------------------
    % PREDICTOR STEP (first-order, identical to original scheme)
    % ------------------------------------------------------------------

    % B predictor
    B_p = (B(n)/h + P.lambda) / (1/h + P.lambda/P.K + P.p*M(n) + P.mu);

    % M predictor
    M_p = (M(n)/h + P.beta * bI_n) / (1/h + P.muM + P.p*B(n));

    % I predictor: boundary condition + PDE along characteristics
    I_p      = zeros(1, ntau);
    I_p(1)   = (1 - CC(n+1)) * P.p * B_p * M_p;
    I_p(2:end) = (I_curr(1:end-1) / h) ./ ...
                 (1/h + P.mu + Gamma2 + sigma_A_n);

    % IG predictor
    IG_p      = zeros(1, ntau);
    IG_p(1)   = CC(n+1) * P.p * B_p * M_p;
    IG_p(2:end) = IG_curr(1:end-1) ./ (1 + h*P.mu + h*GammaG2);

    % G predictor
    G_p = (G(n)/h + bIG_n) / (1/h + P.muG);

    % A predictor
    A_p = (A(n) + h * phi(sumI_n, P.IT, P.s)) / (1 + h * P.muA);

    % ------------------------------------------------------------------
    % CORRECTOR STEP (second-order Crank-Nicolson)
    % ------------------------------------------------------------------

    % Integrals using predicted I_p, IG_p
    bI_p   = h * sum(Gamma1  .* I_p);
    bIG_p  = h * sum(GammaG1 .* IG_p);
    sumI_p = h * sum(I_p);

    % --- B corrector (backward-Euler decay + averaged source, second-order) ---
    % dB/dx = lambda - (lambda/K + mu + p*M)*B
    % Decay term: fully implicit at n+1 (unconditional positivity for any h)
    % Source term: averaged lambda at n and n+1 (both equal lambda, so no change)
    % Nonlinear M coupling: use predictor M_p as estimate of M[n+1]
    % => B[n+1]*(1/h + lambda/K + mu + p*M_p) = B[n]/h + lambda
    % Second-order accuracy: M_p = M[n+1] + O(h), so error in denominator is O(h^2)
    B(n+1) = (B(n)/h + P.lambda) / (1/h + P.lambda/P.K + P.mu + P.p*M_p);

    % --- M corrector (backward-Euler decay + averaged source, second-order) ---
    % dM/dx = beta*bI(x) - (muM + p*B)*M
    % Source beta*bI: averaged using predictor bI_p (trapezoidal, second-order)
    % Decay (muM + p*B)*M: fully implicit at n+1 (positivity for any h, avoids
    %   sign flip in explicit coefficient when h*muM > 2, e.g. h > 0.24 hr)
    % => M[n+1]*(1/h + muM + p*B[n+1]) = M[n]/h + beta/2*(bI_n + bI_p)
    M(n+1) = (M(n)/h + P.beta/2 * (bI_n + bI_p)) / (1/h + P.muM + P.p*B(n+1));

    % --- I corrector: Crank-Nicolson along characteristics ---
    % I[n+1,j+1] = I[n,j] * (2/h - alpha_n_j) / (2/h + alpha_p_j1)
    % alpha = mu + Gamma + sigma*(1-exp(-theta*A))
    sigma_A_p = P.sigma * (1 - exp(-P.theta * A_p));

    alpha_n = P.mu + Gamma1(1:end-1) + sigma_A_n;   % at (n, j),   j=1..ntau-1
    alpha_p = P.mu + Gamma2          + sigma_A_p;   % at pred (n+1, j+1)

    I_next      = zeros(1, ntau);
    I_next(1)   = (1 - CC(n+1)) * P.p * B(n+1) * M(n+1);  % BC uses corrected values
    I_next(2:end) = I_curr(1:end-1) .* (2/h - alpha_n) ./ (2/h + alpha_p);

    % --- IG corrector: Crank-Nicolson along characteristics (no A dependence) ---
    alphaG_n = P.mu + GammaG1(1:end-1);  % at (n, j)
    alphaG_p = P.mu + GammaG2;           % at pred (n+1, j+1)

    IG_next      = zeros(1, ntau);
    IG_next(1)   = CC(n+1) * P.p * B(n+1) * M(n+1);
    IG_next(2:end) = IG_curr(1:end-1) .* (2/h - alphaG_n) ./ (2/h + alphaG_p);

    % --- G corrector (backward-Euler decay + averaged source) ---
    % dG/dx = bIG - muG*G
    % => G[n+1]*(1/h + muG) = G[n]/h + (bIG_n + bIG_p)/2
    G(n+1) = (G(n)/h + (bIG_n + bIG_p)/2) / (1/h + P.muG);

    % --- A corrector (backward-Euler decay + averaged source) ---
    % dA/dx = phi(sumI) - muA*A  (muA = 0 in baseline, so pure integration)
    % => A[n+1]*(1/h + muA) = A[n]/h + (phi_n + phi_p)/2
    A(n+1) = (A(n)/h + (phi(sumI_n, P.IT, P.s) + phi(sumI_p, P.IT, P.s))/2) / (1/h + P.muA);

    % ------------------------------------------------------------------
    % Early exit: stop once G has peaked and dropped back below G_stop.
    % The flag G_peaked prevents an exit during the initial build-up phase
    % when G is still rising from zero.
    % ------------------------------------------------------------------
    if G(n+1) >= G_stop
        G_peaked = true;
    end
    if G_peaked && G(n+1) < G_stop && n >= n_warmup
        break;
    end

    % ------------------------------------------------------------------
    % Advance rolling window
    % ------------------------------------------------------------------
    I_curr  = I_next;
    IG_curr = IG_next;
    infection_level(n+1) = h * sum(I_curr);

end
end
