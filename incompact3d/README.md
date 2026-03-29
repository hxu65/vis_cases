# Numerical Dissipation in Under-Resolved High-Order Simulations of the Taylor-Green Vortex

## Introduction

High-order numerical methods are widely used in computational fluid dynamics for their superior accuracy at resolved scales. However, this accuracy comes with a trade-off: high-order compact finite-difference schemes, such as those employed in Xcompact3d, are designed to minimize numerical dissipation. While this property is desirable for well-resolved simulations, it becomes problematic when the computational grid is too coarse to capture the smallest dynamically active scales of turbulence.

In under-resolved configurations, the lack of numerical dissipation allows unresolved kinetic energy to accumulate near the grid cutoff wavenumber — a phenomenon known as spectral energy pile-up. This pile-up manifests as spurious velocity overshoots, growing pre-projection divergence, and ultimately degraded or non-physical solutions. Understanding this mechanism is essential for practitioners who must choose between fully resolved DNS, explicit subgrid-scale (SGS) modelling, and implicit LES strategies.

This case study reproduces the under-resolution problem using the three-dimensional Taylor-Green vortex (TGV) at Reynolds numbers of 5,000 and 10,000 on a deliberately coarse 65^3 grid, demonstrating how spectral energy pile-up develops when the grid spacing far exceeds the Kolmogorov scale.

## Background

### High-Order Compact Schemes and Numerical Dissipation

Xcompact3d solves the incompressible Navier-Stokes equations using 6th-order compact finite-difference schemes for spatial differentiation. These schemes achieve spectral-like resolution — the modified wavenumber closely follows the exact differentiation curve across most of the resolvable wavenumber range. Crucially, the schemes are purely centred, introducing no intrinsic upwind dissipation. This means that energy at poorly resolved wavenumbers near the grid cutoff is neither damped nor removed by the numerical method itself.

For DNS on adequately refined grids (where the mesh spacing dx is comparable to the Kolmogorov scale eta), the energy spectrum at the cutoff wavenumber is already exponentially small, and the absence of numerical dissipation poses no difficulty. For under-resolved computations, however, the energy at the cutoff remains significant, and the truncation error acts to redistribute energy into the highest resolvable modes rather than dissipating it.

### The Implicit LES Approach

Dairay et al. (2017) proposed an alternative to explicit SGS models: introducing a controlled, targeted numerical dissipation through the discretization of the viscous term. This approach, termed *controlled implicit LES*, modifies the effective viscosity using a spectral vanishing viscosity (SVV) kernel parameterized by:

- **nu0nu** (`nu0/nu`): the ratio of the maximum artificial viscosity to the molecular viscosity, controlling the amplitude of the added dissipation.
- **cnu**: a parameter governing the spectral shape of the SVV kernel, determining at which wavenumbers the artificial dissipation acts.

The key insight is that by calibrating these parameters using a Pao-like spectral closure, the numerical dissipation can be made physically consistent — concentrating the added dissipation near the cutoff wavenumber where it is needed, while preserving the accuracy of the scheme at well-resolved scales.

### Reference

> Dairay, T., Lamballais, E., Laizet, S., & Vassilicos, J. C. (2017). Numerical dissipation vs. subgrid-scale modelling for large eddy simulation. *Journal of Computational Physics*, 337, 252–274.

The paper demonstrates that this implicit approach outperforms both standard and dynamic Smagorinsky models in terms of agreement with DNS, while also achieving numerical convergence — a property that explicit SGS models on the same grid cannot guarantee. A copy of the paper is included as `2017_laizet_jcp.pdf`.

## Methodology

### Test Case: Taylor-Green Vortex

The TGV is a canonical benchmark in which a smooth, deterministic initial condition

```
u(x, y, z, 0) =  sin(x) cos(y) cos(z)
v(x, y, z, 0) = -cos(x) sin(y) cos(z)
w(x, y, z, 0) =  0
```

evolves through vortex stretching and instability into fully developed turbulence within a triply periodic domain of size (2pi)^3. The flow exhibits a well-characterized enstrophy peak at t ~ 8-10 driven by vortex stretching, followed by viscous decay. Because the evolution is entirely deterministic (no random forcing or boundary effects), discrepancies from reference solutions can be unambiguously attributed to numerical errors.

### Simulation Parameters

Both cases use identical numerical settings except for the Reynolds number:

| Parameter | Scenario 1 | Scenario 2 |
|-----------|-----------|-----------|
| **Directory** | `scenario1_tgv_underresolved/` | `tgv_Re10000_65/` |
| **Re** | 5,000 | 10,000 |
| **Grid** | 65 x 65 x 65 | 65 x 65 x 65 |
| **dx** | 0.0491 | 0.0491 |
| **Kolmogorov scale (eta)** | ~0.0053 | ~0.0028 |
| **dx / eta** | ~9x | ~17x |
| **dt** | 0.001 | 0.001 |
| **Time integration** | RK3 | RK3 |
| **Spatial scheme** | 6th-order compact | 6th-order compact |
| **nu0nu** | 4.0 (standard) | 4.0 (standard) |
| **cnu** | 0.44 (standard) | 0.44 (standard) |
| **Total steps** | 8,000 (t = 8.0) | 10,000 (t = 10.0) |
| **ilesmod** | 0 (off) | 0 (off) |

In both configurations, the implicit LES module is disabled (`ilesmod=0`), and `nu0nu = 4.0` with `cnu = 0.44` correspond to the standard 6th-order scheme without enhanced dissipation. For comparison, the ILES configurations studied by Dairay et al. (2017) used `nu0nu = 63` (steep SVV) or `nu0nu = 115` (sharp SVV) at Re = 5,000 on a 160^3 grid.

## Results

### Scenario 1: Re = 5,000 on 65^3

The simulation completed all 8,000 time steps (t = 8.0) in approximately 69 minutes of wall-clock time. However, the solution shows clear signatures of under-resolution:

| Time (t) | DIV U\* max | U max | W max | CFL_z |
|---------:|----------:|------:|------:|------:|
| 0.001 | 6.65e-04 | 1.000 | 0.000 | 0.000 |
| 1.0 | 3.35e-04 | 0.961 | 0.260 | 0.005 |
| 2.0 | 3.68e-04 | 0.936 | 0.559 | 0.011 |
| 3.0 | 1.71e-03 | 0.943 | 0.841 | 0.017 |
| 4.0 | 2.77e-02 | 0.901 | 0.951 | 0.019 |
| 4.5 | 1.24e-01 | 1.160 | 1.240 | 0.025 |
| 5.0 | 1.22e-01 | 1.095 | 1.170 | 0.024 |
| 6.0 | 2.08e-01 | 1.122 | 1.111 | 0.023 |
| 7.0 | 3.42e-01 | 1.751 | 1.227 | 0.025 |
| 8.0 | 3.99e-01 | 1.385 | 1.496 | 0.031 |

The pre-projection divergence (DIV U\*) — which measures how far the intermediate velocity field deviates from the divergence-free constraint before the pressure correction — grows by nearly three orders of magnitude, from O(10^-4) at t = 0 to a peak of **0.613 at t = 6.8**. Velocity magnitudes exceed the initial maximum of 1.0 starting at t ~ 4.5, with the peak reaching **U_max = 1.785** at t = 6.9 (a 79% overshoot). These spurious velocities are a direct consequence of energy piling up at the grid scale.

The post-projection divergence (DIV U) remains at machine precision (~10^-13) throughout, confirming that the Poisson solver functions correctly — the problem lies in the quality of the velocity field itself.

### Scenario 2: Re = 10,000 on 65^3

At Re = 10,000, the Kolmogorov scale shrinks to eta ~ 0.0028, making the grid-to-dissipation-scale ratio approximately 17:1. The simulation reached step 8,575 (t = 8.575) before being terminated by a SLURM node failure (not numerical blow-up).

| Time (t) | DIV U\* max | U max | W max | CFL_z |
|---------:|----------:|------:|------:|------:|
| 0.001 | 6.65e-04 | 1.000 | 0.000 | 0.000 |
| 1.0 | 3.35e-04 | 0.962 | 0.260 | 0.005 |
| 2.0 | 3.70e-04 | 0.937 | 0.560 | 0.011 |
| 3.0 | 1.72e-03 | 0.946 | 0.847 | 0.017 |
| 4.0 | 3.20e-02 | 0.913 | 0.963 | 0.020 |
| 4.5 | 1.32e-01 | 1.188 | 1.256 | 0.026 |
| 5.0 | 1.49e-01 | 1.124 | 1.181 | 0.024 |
| 6.0 | 3.06e-01 | 1.314 | 1.145 | 0.023 |
| 7.0 | 3.66e-01 | 1.963 | 1.437 | 0.029 |
| 8.0 | 6.20e-01 | 1.478 | 1.579 | 0.032 |
| 8.5 | 5.97e-01 | 1.315 | 1.392 | 0.028 |

The degradation follows the same qualitative pattern but is quantitatively more severe:

| Metric | Re = 5,000 | Re = 10,000 | Ratio |
|--------|----------:|----------:|------:|
| Peak DIV U\* max | 0.613 (t = 6.8) | 0.945 (t = 8.1) | 1.54x |
| Peak U max | 1.785 (t = 6.9) | 2.173 (t = 6.9) | 1.22x |
| Peak W max | 1.691 (t = 7.6) | 1.879 (t = 7.4) | 1.11x |

### Key Observations

1. **Onset timing is Re-independent.** Both cases show the onset of significant degradation at t ~ 4.0-4.5, consistent with the TGV transition time being only weakly dependent on Reynolds number. The early-time evolution (t < 3) is virtually identical between the two cases, as the flow remains laminar and well-resolved at this grid spacing.

2. **Severity scales with under-resolution.** The peak pre-projection divergence at Re = 10,000 is 54% larger than at Re = 5,000, and spurious velocity overshoots reach 117% above the initial maximum (compared to 79% at Re = 5,000). This is consistent with the expectation that a larger dx/eta ratio produces more severe aliasing.

3. **The simulations survive but produce non-physical solutions.** Neither case diverges numerically — CFL numbers remain well below unity (peak ~ 0.04), and the pressure projection maintains a divergence-free field at machine precision. The problem is not blow-up but degraded accuracy: the computed velocity field contains energy at scales the grid cannot represent, producing spurious small-scale fluctuations that contaminate the large-scale solution.

4. **Comparison with Dairay et al. (2017).** In Table 1 of the reference paper, Case 2 (no-model LES at Re = 5,000 on a 160^3 grid, dx = 0.039) already shows an "unrealistic pile-up of energy in the vicinity of the cutoff wavenumber." The present cases use a 65^3 grid (dx = 0.049), which is even coarser, producing a more pronounced version of the same phenomenon.

## Remediation Strategies

Two approaches can address the under-resolution problem:

1. **Increase grid resolution (DNS approach).** Refining the grid to 128^3, 257^3, or finer brings dx closer to eta, reducing the aliasing error. Dairay et al. used 1280^3 for DNS at Re = 5,000 and 2048^3 at Re = 10,000. This approach is computationally expensive but introduces no modelling error.

2. **Enable controlled implicit LES (ILES approach).** By increasing `nu0nu` from 4.0 to values in the range of 47-115 (depending on the SVV kernel shape), a targeted dissipation is added near the cutoff wavenumber. Combined with appropriate `cnu` values, this removes the spectral pile-up while preserving accuracy at well-resolved scales. Dairay et al. (2017) showed that this approach on a 160^3 grid at Re = 5,000 produces closer agreement with DNS than either the standard or dynamic Smagorinsky models.

## File Structure

```
incompact3d/
├── README.md                          # This document
├── 2017_laizet_jcp.pdf                # Reference paper (Dairay et al. 2017)
├── images/                            # Simulation visualization recording
├── scenario1_tgv_underresolved/       # Re = 5,000, 65^3, no ILES
│   ├── input.i3d                      # Xcompact3d configuration
│   ├── output.log                     # Full simulation output
│   ├── decomp_2d_setup.log            # Domain decomposition log
│   └── restart.info                   # Checkpoint metadata
└── tgv_Re10000_65/                    # Re = 10,000, 65^3, no ILES
    ├── input.i3d                      # Xcompact3d configuration
    ├── output.log                     # Partial output (node failure at t=8.575)
    └── decomp_2d_setup.log            # Domain decomposition log
```

## Video

A visualization recording of the simulation is available in [`images/`](images/).
