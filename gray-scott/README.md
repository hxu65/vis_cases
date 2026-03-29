# Gray-Scott Reaction-Diffusion: Chaotic Divergence under MPI Decomposition

Visualization cases demonstrating how the Gray-Scott reaction-diffusion system in its spatiotemporal chaos regime produces qualitatively different pattern states when run with different MPI decompositions — despite identical parameters and initial conditions.

## Background

The Gray-Scott model simulates two chemical species (U and V) reacting and diffusing on a 3D grid:

```
∂U/∂t = Du ∇²U − UV² + F(1 − U)
∂V/∂t = Dv ∇²V + UV² − (F + k)V
```

With parameters in the **spatiotemporal chaos regime** (F=0.02, k=0.048), the system has no stable attractor. Tiny floating-point differences — such as those introduced by changing the order of MPI reduction operations across different domain decompositions — are amplified chaotically over time. After enough simulated steps, two runs with different process layouts diverge into completely different pattern states (e.g., spots vs. labyrinths), with no ground truth to say which is "correct."

## Source Code

The simulation uses the ADIOS2 Gray-Scott example, which ships with every ADIOS2 release:

- **C++ (ADIOS2):** [github.com/ornladios/ADIOS2](https://github.com/ornladios/ADIOS2) (`examples/` directory)
- **Standalone:** [github.com/keichi/gray-scott](https://github.com/keichi/gray-scott)
- **Julia port:** [github.com/JuliaORNL/GrayScott.jl](https://github.com/JuliaORNL/GrayScott.jl) — supports Float32 vs Float64 precision and CPU/CUDA/AMDGPU backends via the same settings JSON

Configs are pinned to ADIOS2 v2.10.x tagged releases. The Julia version is archived on Zenodo.

## Experiment Design

### Core idea

Run the **identical** chaos-regime configuration with different MPI process layouts. The floating-point reduction order changes with the decomposition, seeding divergence that the chaotic dynamics amplify into macroscopically different patterns.

### Default parameters (chaos regime)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `F`       | 0.02  | Feed rate |
| `k`       | 0.048 | Kill rate |
| `Du`      | 0.2   | Diffusion rate of U |
| `Dv`      | 0.1   | Diffusion rate of V |
| `dt`      | 1.0   | Time step |
| `L`       | 64    | Grid size (L³) |
| `noise`   | 0.01  | Initial perturbation amplitude |

### MPI decomposition runs

| Run | Processes | Frames | Description |
|-----|-----------|--------|-------------|
| `run_2x1x1` | 2  | 80  | 2 processes along X |
| `run_2x2x1` | 4  | 80  | 2x2 slab in XY |
| `run_2x2x2` | 8  | 100 | Full 3D decomposition |
| `run_4x2x1` | 8  | 100 | Asymmetric 4x2 slab in XY |

All runs use `settings-chaos.json` (F=0.02, k=0.048). The initial conditions and parameters are identical — only the MPI topology differs.

### What to look for

- **Early time steps:** All decompositions show nearly identical volume renders. MPI partition boundaries are visible as seam lines on the isosurface.
- **Late time steps:** The `2x2x2` and `4x2x1` runs (8 processes) develop dramatically fragmented, chaotic isosurface structures with mixed spots and labyrinthine features. The `2x1x1` and `2x2x1` runs (2–4 processes) evolve into qualitatively different, more coherent shapes. No two decompositions converge to the same final state.

This is the key result: **there is no "correct" output** — the chaos regime lacks a stable attractor, so any pattern is equally valid.

## Configuration Files

| Config | F | k | Du | dt | L | Purpose |
|--------|---|---|----|----|---|---------|
| `settings-chaos.json` | 0.02 | 0.048 | 0.2 | 1.0 | 64 | **Primary experiment** — spatiotemporal chaos regime |
| `settings-analysis.json` | 0.02 | 0.048 | 0.2 | 1.0 | 64 | Same chaos regime, longer run (5000 steps) for analysis |
| `settings-stable.json` | 0.04 | 0.06 | 0.2 | 1.0 | 64 | Stable spots regime — control case |
| `settings-unstable-params.json` | 0.08 | 0.03 | 0.2 | 1.0 | 64 | Different F/k — tests parameter sensitivity |
| `settings-unstable-Du.json` | 0.04 | 0.06 | 0.3 | 1.0 | 64 | Elevated Du (0.3) — tests diffusion sensitivity |
| `settings-unstable-dt.json` | 0.04 | 0.06 | 0.2 | 2.0 | 64 | Doubled time step (dt=2.0) — tests temporal stability |
| `settings-files.json` | 0.01 | 0.05 | 0.2 | 2.0 | 128 | Larger grid (128³), different regime |
| `settings-files2.json` | 0.01 | 0.05 | 0.2 | 2.0 | 128 | Same as above, alternate output path (`gs2.bp`) |

### Stability comparison

The `settings-stable.json` config (F=0.04, k=0.06) sits in a regime with a stable attractor — different MPI decompositions should converge to the same pattern. This serves as a **control** to confirm that divergence in the chaos runs is physical (chaotic sensitivity), not a software bug.

The three `unstable-*` configs perturb individual parameters away from the stable baseline to explore where instability onset begins.

## Reproducing

```bash
# Build ADIOS2 with the Gray-Scott example enabled
# Then run with different MPI layouts:
mpirun -np 2 gray-scott settings-chaos.json   # 2x1x1
mpirun -np 4 gray-scott settings-chaos.json   # 2x2x1
mpirun -np 8 gray-scott settings-chaos.json   # 2x2x2 or 4x2x1
```

The process layout is determined by the MPI decomposition logic in the code. To force a specific layout (e.g., 4x2x1 vs 2x2x2 with 8 processes), set the decomposition explicitly via command-line arguments — see the ADIOS2 Gray-Scott README for details.

## Additional deviation experiment

Beyond MPI decomposition, the Julia port enables a **precision deviation** experiment: run the same F=0.02, k=0.048 config with Float32 vs Float64. The reduced precision introduces rounding differences that, like MPI reordering, are chaotically amplified into completely different final states — another demonstration that no single "correct" output exists in this regime.
