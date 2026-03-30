# withinhost_higher_order

Second-order MATLAB solver for the age-structured within-host malaria PDE model from Patterson et al. (2026) *Evolution*.

Companion repository to [malaria_within_host](https://github.com/patterd2/malaria_within_host), providing a faster, more accurate numerical scheme.

## What's new

The original solver uses a first-order implicit upwind scheme (O(h) accuracy), requiring h = 0.125 hr for sufficient accuracy. This repo implements:

1. **Second-order predictor-corrector** — Crank-Nicolson along characteristics for the age-structured PDEs (I, IG) and semi-implicit trapezoidal updates for the ODEs (B, M, G, A). Achieving O(h²) accuracy allows h = 0.25 hr with comparable accuracy to the original h = 0.125 hr (2× fewer time steps).

2. **Rolling-window memory** — The full I(nx×ntau) and IG(nx×ntau) matrices (~12 GB at baseline parameters) are replaced by a pair of O(ntau) row vectors. This eliminates cache thrashing and drastically reduces memory usage.

Combined, these changes give roughly **4–10× speedup** for single simulations and optimization runs.

## Files

| File | Description |
|------|-------------|
| `within_host_model_2nd_order.m` | Core solver (2nd-order predictor-corrector) |
| `run_full_model.m` | Run a single simulation |
| `withinhost_model_optimization.m` | Fitness objective for strategy optimization |
| `run_strategy_optimization.m` | Optimize plastic investment strategy via Nelder-Mead |
| `convergence_test.m` | Verify O(h²) convergence vs O(h) for original scheme |
| `standard_plotting_2nd_order.m` | Plotting script (called by run_full_model.m) |
| `baseline_parameter_set.m` | Model parameters (identical to original repo) |
| `gamma_fun.m` | Burst hazard rate for infected RBCs |
| `gamma_G.m` | Maturation hazard rate for developing gametocytes |
| `phi.m` | Immune activation sigmoid function |
| `betaHV.m` | Host-to-vector infectivity function |
| `simps.m` | Simpson's rule numerical integration |

## Usage

### Single simulation

```matlab
run_full_model
```

Set `h = 0.25` (instead of `0.125`) in `run_full_model.m` for 2× speedup.

### Strategy optimization

```matlab
run_strategy_optimization
```

### Convergence verification

```matlab
convergence_test
```

Produces a log-log plot showing slope ≈ 1 (original) and ≈ 2 (new scheme).

## Numerical scheme summary

For each step n → n+1 with step size h (in both x and τ):

**Predictor** (first-order, same as original):
```
B*, M*, I*, IG*, G*, A*  ←  current values at step n
```

**Corrector** (Crank-Nicolson, second-order):
```
B[n+1] = (B[n]·(2/h − λ/K − μ) + 2λ − p·(B[n]·M[n] + B*·M*)) / (2/h + λ/K + μ)

I[n+1,j+1] = I[n,j] · (2/h − α_n[j]) / (2/h + α_p[j+1])
  where α = μ + Γ(τ) + σ·(1 − exp(−θA))   (positivity guaranteed since α_max ≈ 1.02 ≪ 2/h)

A[n+1] = (A[n]·(2/h − μA) + φ(∫I_n dτ) + φ(∫I* dτ)) / (2/h + μA)
```

## Reference

Patterson DD, Childs LM, Cao N, Srivathsa-Abber S, Greischar MA (2026).
*Immunity can impose a reproduction-survival tradeoff on human malaria parasites.*
Evolution 80(2): 396–411. https://doi.org/10.1093/evolut/qpaf238
