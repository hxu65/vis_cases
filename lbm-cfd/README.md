# LBM-CFD: Lattice-Boltzmann Method 2D Computational Fluid Dynamics

A parallel 2D Lattice-Boltzmann CFD solver with integrated in-situ visualization, interactive web-based steering, and offline post-processing capabilities.

## Overview

This project implements a **D2Q9 Lattice-Boltzmann Method** (LBM) simulation for modeling viscous fluid flow through a channel with configurable obstacles. The simulation is parallelized using **MPI** with 2D domain decomposition and features a three-layer visualization architecture coupling [Ascent](https://ascent.readthedocs.io/) for in-situ rendering, [Trame](https://kitware.github.io/trame/) for interactive web visualization, and [ParaView](https://www.paraview.org/) for offline analysis.

### Physical Setup

The default configuration models **corn syrup (25 C) flowing at 0.75 m/s through a 2 m pipe** over 8 seconds of physical time:

| Parameter | Value |
|---|---|
| Fluid density | 1380 kg/m^3 |
| Dynamic viscosity | 1.3806 Pa s |
| Inlet velocity | 0.75 m/s |
| Domain resolution | 600 x 240 cells |
| Simulation timesteps | 20,000 (default) |

Two pairs of vertical barriers with center gaps create flow restrictions, producing vortex shedding patterns characteristic of low-Reynolds-number flows.

## Simulation Results

Vorticity field evolution showing vortex formation at barrier gaps, rendered with a divergent (blue-red) colormap:

| Early flow development | Vortex evolution | Onset of instability |
|---|---|---|
| ![Step 10](images/step%2010.png) | ![Step 15](images/step%2015.png) | ![Step 20](images/step%2020.png) |

## Method

### Lattice-Boltzmann D2Q9 Scheme

The solver uses the BGK (Bhatnagar-Gross-Krook) single-relaxation-time collision operator with 9 discrete velocity directions (1 rest + 4 cardinal + 4 diagonal). The relaxation parameter is computed as:

```
omega = 1 / (3 * viscosity + 0.5)
```

Equilibrium distribution weights follow the standard D2Q9 model: 4/9 (rest), 1/9 (cardinal), 1/36 (diagonal).

Each timestep consists of:
1. **Collision** -- Relax distributions toward local equilibrium (BGK approximation)
2. **Streaming** -- Propagate distributions to neighboring lattice sites
3. **Bounce-back** -- Reflect distributions at solid boundaries (barrier obstacles)

### MPI Parallelization

The domain is decomposed into a 2D grid of MPI ranks using an optimal factorization strategy. Each rank maintains one-cell overlap (halo) regions exchanged via custom MPI datatypes for efficient column and row transfers. All 14 scalar fields per cell (9 distribution functions, density, 2 velocity components, vorticity, speed) are allocated as a single contiguous array for cache-friendly access.

### Stability Monitoring and Recovery

The simulation includes an automatic stability detection and checkpoint/rollback mechanism:
- Density fields are monitored at regular intervals via `MPI_Allreduce`
- State checkpoints are saved every 500 steps
- On instability detection, the simulation reverts to the last checkpoint and doubles the timestep count (up to 5 recovery attempts)
- A `--force-unstable` mode is available for testing with intentionally reduced timesteps

## Visualization Architecture

```
Simulation (C++/MPI)
    |
    v  Ascent in-situ
    |
    +--[stable]--> ascent_trame_bridge.py --> Trame Web UI
    |                                           |
    |              steering (barriers, speed) <--+
    |
    +--[unstable]--> ascent_rescue.py --> email alert + checkpoint rollback
    |
    +--[--output-vts]--> VTS files --> create_pvd_from_vts.py --> ParaView
```

### In-Situ Visualization (Ascent)

When built with Ascent support, the simulation publishes vorticity fields as Conduit Blueprint meshes at configurable intervals. Ascent executes Python extract scripts that route data to either the interactive Trame path (stable state) or the rescue/notification path (unstable state).

### Interactive Web Visualization (Trame)

The Trame application provides a real-time web interface with:
- **Live vorticity rendering** with selectable colormaps (divergent, turbo, inferno)
- **Flow speed control** via slider (0.25 -- 1.50 m/s)
- **Interactive barrier editing** -- click-and-drag to place vertical/horizontal obstacles
- **Computational steering** -- parameter changes propagate back to the running simulation via a queue-based callback architecture

### Offline Post-Processing (ParaView)

Running with `--output-vts` exports VTK StructuredGrid files per timestep containing vorticity (clamped and raw) and barrier mask fields. The `create_pvd_from_vts.py` utility bundles these into PVD time-series collections for ParaView animation.

## Build

**With Ascent support:**
```bash
env ASCENT_DIR=<ascent_install_dir>/install/ascent-checkout make
```

**Without Ascent support (VTS output only):**
```bash
make
```

Requires `mpicxx` with C++14 support. Cross-platform build system supports Linux, macOS, and Windows (via MSVC/MSMPI).

## Usage

### Run the Trame Server

Activate the Python virtual environment created for the Ascent installation, then:

```bash
python trame/trame_app.py --host 0.0.0.0 --port <port> --server --timeout 0
```

### Run the Simulation

Ensure the Trame server is running and the web page is open, then:

```bash
PYTHON_SITE_PKG="<venv_path>/lib/python3.12/site-packages"
ASCENT_DIR="<ascent_install_dir>/install"
export PYTHONPATH=$PYTHONPATH:$PYTHON_SITE_PKG:$ASCENT_DIR/ascent-checkout/python-modules/:$ASCENT_DIR/conduit-v0.9.2/python-modules/

mpiexec -np <num_procs> ./bin/lbmcfd [options]
```

**Command-line options:**

| Flag | Description |
|---|---|
| `--steps <N>` | Override total simulation timesteps (default: 20,000) |
| `--output-vts` | Export VTS files to disk for ParaView |
| `--output-dir <path>` | Directory for VTS output (default: `paraview/`) |
| `--force-unstable` | Run with reduced timesteps to trigger instability |

### Generate ParaView Collections

```bash
make pvd
# or manually:
python create_pvd_from_vts.py --input-dir paraview
```

## Project Structure

```
lbm-cfd/
├── src/main.cpp              # Simulation loop, Ascent/Trame integration
├── include/lbmd2q9_mpi.hpp   # LBM D2Q9 solver with MPI parallelization
├── ascent/
│   ├── ascent_trame_bridge.py  # In-situ to Trame data bridge (stable path)
│   └── ascent_rescue.py        # Instability handler with email alerts
├── trame/trame_app.py        # Interactive web visualization server
├── resrc/                    # Colormap PNG resources (divergent, turbo, inferno)
├── create_pvd_from_vts.py    # VTS-to-PVD collection generator
├── Makefile                  # Multi-platform build system
└── images/                   # Visualization snapshots
```

## Dependencies

- **C++ compiler** with C++14 and MPI support (`mpicxx`)
- **Ascent** (optional) -- in-situ visualization framework with Conduit
- **Python 3.12+** with packages: `trame`, `numpy`, `opencv-python`
- **ParaView** (optional) -- for offline VTS/PVD analysis
