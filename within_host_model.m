function [B, M, I, IG, G, A, CC] = within_host_model(h, t0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC)

global P

tau = (0:h:tau_max)';
ntau = length(tau);
x = (t0:h:X_max)';
nx = length(x); % all discretization parameters chosen equal

% allocation of memory
B = NaN(nx,1); % B(x)
M = NaN(nx,1); % M(x)
I = NaN(nx,ntau); % I(x,tau)
IG = NaN(nx,ntau); % IG(x,tau)
G = NaN(nx,1); % G(x)
A = NaN(nx,1); % A(x)

% assisgn initial conditions
B(1) = B0;
M(1) = M0;
I(1,:) = I0; % x = 0, I(0,tau)
IG(1,:) = IG0; % x = 0, IG(x,tau)
G(1) = G0;
A(1) = A0;
%I(1,1) = (1-P.c)*P.p*B0*M0;
%IG(1,1) = P.c*P.p*B0*M0;
%% time since infection evolution
% precalculate some vectors to avoid many calls to gamcdf
Gamma1 = gamma_fun(tau,h);
Gamma2 = gamma_fun(tau(2:end)',h);
Gamma3 = h*gamma_G(tau(2:end),h)';
Gamma4 = gamma_G(tau,h);
for n = 1:nx-1 % evolving on the time since infection time scale
    % evolve ODEs for red blood cells and merozoites
    %fprintf('%i \n',n);
    B(n+1) = (B(n)/h + P.lambda)./(1/h + P.lambda/P.K + P.p*M(n) + P.mu);
    M(n+1) = (M(n)/h + h*P.beta*sum(Gamma1.*I(n,:)'))./( 1/h + P.muM +P.p*B(n));
    % BC for I
    I(n+1,1) = (1-CC(n+1))*P.p*B(n+1)*M(n+1);
    % evolve I
    I(n+1,2:end) = (I(n,1:end-1)/h)./(1/h + P.mu + Gamma2 + P.sigma.*(1-exp(-P.theta*A(n))) );
    % BC for IG
    IG(n+1,1) = CC(n+1)*P.p*B(n+1)*M(n+1);
    % evolve IG
    IG(n+1,2:end) = IG(n,1:end-1)./(1 + h*P.mu + Gamma3);
    % evolve ODEs for gametocytes and the immune activation level
    G(n+1) = ((1/h)*G(n) + h*sum(Gamma4.*IG(n,:)'))./(1/h + P.muG);
    %A(n+1) = (A(n) + h*(phi( h*sum(I(n,:),2), P.IT, P.s)))/(1 + h*P.muA); 
    A(n+1) = (A(n) + h*phi(h*sum(I(n,:),2),P.IT,P.s))/(1 + h*P.muA);
end
end
