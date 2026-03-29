#!/usr/bin/env python3
"""XY slice comparison for dense LJ fluid: correct vs wrong timestep.
Uses local density, coordination number, and Voronoi analysis to make
structural differences VISIBLE in the slice."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial.distance import cdist

BASE = "/mnt/common/hxu40/File/lammps_cases/case1_timestep_explosion"

def read_xyz_frame(filename, frame_idx):
    with open(filename) as f:
        lines = f.readlines()
    i, current = 0, 0
    last_start = 0
    while i < len(lines):
        natoms = int(lines[i].strip())
        if frame_idx >= 0 and current == frame_idx:
            coords = []
            for j in range(i + 2, i + 2 + natoms):
                parts = lines[j].split()
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            return np.array(coords)
        last_start = i
        i += natoms + 2
        current += 1
    # last frame
    natoms = int(lines[last_start].strip())
    coords = []
    for j in range(last_start + 2, last_start + 2 + natoms):
        parts = lines[j].split()
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.array(coords)

def count_frames(filename):
    with open(filename) as f:
        lines = f.readlines()
    i, count = 0, 0
    while i < len(lines):
        try:
            natoms = int(lines[i].strip())
            i += natoms + 2
            count += 1
        except:
            break
    return count

def xy_slice(coords, z_frac=0.08):
    zmid = np.median(coords[:, 2])
    zrange = coords[:, 2].max() - coords[:, 2].min()
    z_lo = zmid - z_frac * zrange
    z_hi = zmid + z_frac * zrange
    mask = (coords[:, 2] >= z_lo) & (coords[:, 2] <= z_hi)
    return coords[mask]

def coordination_number(coords_2d, r_cut=1.5):
    """Count neighbors within r_cut for each atom in 2D."""
    d = cdist(coords_2d[:, :2], coords_2d[:, :2])
    np.fill_diagonal(d, np.inf)
    return np.sum(d < r_cut, axis=1)

def local_density(coords_2d, r_probe=2.0):
    """Local 2D density around each atom."""
    d = cdist(coords_2d[:, :2], coords_2d[:, :2])
    np.fill_diagonal(d, np.inf)
    count = np.sum(d < r_probe, axis=1)
    return count / (np.pi * r_probe**2)

def nearest_neighbor_dist(coords):
    d = cdist(coords, coords)
    np.fill_diagonal(d, np.inf)
    return np.min(d, axis=1)

# Load data
nf_c = count_frames(f"{BASE}/trajectory_dense_correct.xyz")
nf_w = count_frames(f"{BASE}/trajectory_dense_wrong.xyz")
print(f"Correct: {nf_c} frames, Wrong: {nf_w} frames")

# ============================================================
# Plot 1: XY slice colored by COORDINATION NUMBER
# ============================================================
print("Generating coordination number comparison...")

fig, axes = plt.subplots(2, 4, figsize=(24, 12))

frame_pcts = [0, 0.33, 0.66, 1.0]
labels = ['Start', '33%', '66%', 'End']

for col, (pct, label) in enumerate(zip(frame_pcts, labels)):
    idx_c = min(int(pct * (nf_c - 1)), nf_c - 1)
    idx_w = min(int(pct * (nf_w - 1)), nf_w - 1)

    coords_c = read_xyz_frame(f"{BASE}/trajectory_dense_correct.xyz", idx_c)
    coords_w = read_xyz_frame(f"{BASE}/trajectory_dense_wrong.xyz", idx_w)

    sl_c = xy_slice(coords_c, z_frac=0.06)
    sl_w = xy_slice(coords_w, z_frac=0.06)

    cn_c = coordination_number(sl_c, r_cut=1.5)
    cn_w = coordination_number(sl_w, r_cut=1.5)

    # Correct (top)
    ax = axes[0, col]
    sc = ax.scatter(sl_c[:, 0], sl_c[:, 1], s=120, c=cn_c, cmap='coolwarm',
                   vmin=0, vmax=10, alpha=0.8, edgecolors='k', linewidths=0.3)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2)
    ax.set_title(f'{label} (frame {idx_c})\nMean CN={cn_c.mean():.1f}', fontsize=11)
    if col == 0:
        ax.set_ylabel('CORRECT dt=0.005\nY (σ)', fontsize=12, color='blue', fontweight='bold')

    # Wrong (bottom)
    ax = axes[1, col]
    sc = ax.scatter(sl_w[:, 0], sl_w[:, 1], s=120, c=cn_w, cmap='coolwarm',
                   vmin=0, vmax=10, alpha=0.8, edgecolors='k', linewidths=0.3)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2)
    ax.set_title(f'{label} (frame {idx_w})\nMean CN={cn_w.mean():.1f}', fontsize=11)
    ax.set_xlabel('X (σ)', fontsize=11)
    if col == 0:
        ax.set_ylabel('WRONG dt=0.02\nY (σ)', fontsize=12, color='red', fontweight='bold')

# Colorbar
cbar = fig.colorbar(sc, ax=axes, shrink=0.6, label='Coordination Number (r < 1.5σ)')

fig.suptitle('XY Slice: Atoms Colored by Coordination Number (Dense LJ Fluid, ρ*=0.8442)\n'
             'Wrong timestep → heating → lower coordination (more blue = fewer neighbors)',
             fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/dense_coordination_comparison.png", dpi=150)
plt.close()
print("  -> dense_coordination_comparison.png")

# ============================================================
# Plot 2: Side-by-side FINAL frame with multiple diagnostics
# ============================================================
print("Generating final-frame diagnostic panel...")

coords_c = read_xyz_frame(f"{BASE}/trajectory_dense_correct.xyz", -1)
coords_w = read_xyz_frame(f"{BASE}/trajectory_dense_wrong.xyz", -1)

sl_c = xy_slice(coords_c, z_frac=0.08)
sl_w = xy_slice(coords_w, z_frac=0.08)

fig, axes = plt.subplots(2, 3, figsize=(21, 14))

# Row 1: Correct
# Col 0: Raw positions
ax = axes[0, 0]
ax.scatter(sl_c[:, 0], sl_c[:, 1], s=100, c='steelblue', alpha=0.7,
          edgecolors='darkblue', linewidths=0.5)
ax.set_aspect('equal'); ax.grid(True, alpha=0.2)
ax.set_title(f'Raw Positions ({len(sl_c)} atoms)', fontsize=12)
ax.set_ylabel('CORRECT dt=0.005\nY (σ)', fontsize=13, color='blue', fontweight='bold')

# Col 1: Colored by coordination number
cn_c = coordination_number(sl_c, r_cut=1.5)
ax = axes[0, 1]
sc1 = ax.scatter(sl_c[:, 0], sl_c[:, 1], s=100, c=cn_c, cmap='RdYlGn',
                vmin=1, vmax=8, alpha=0.8, edgecolors='k', linewidths=0.3)
ax.set_aspect('equal'); ax.grid(True, alpha=0.2)
ax.set_title(f'Coordination Number\nMean={cn_c.mean():.2f}, Std={cn_c.std():.2f}', fontsize=12)
plt.colorbar(sc1, ax=ax, shrink=0.7, label='CN')

# Col 2: Colored by local density
ld_c = local_density(sl_c, r_probe=1.8)
ax = axes[0, 2]
sc2 = ax.scatter(sl_c[:, 0], sl_c[:, 1], s=100, c=ld_c, cmap='viridis',
                vmin=0.3, vmax=1.5, alpha=0.8, edgecolors='k', linewidths=0.3)
ax.set_aspect('equal'); ax.grid(True, alpha=0.2)
ax.set_title(f'Local Density (r<1.8σ)\nMean={ld_c.mean():.3f}', fontsize=12)
plt.colorbar(sc2, ax=ax, shrink=0.7, label='ρ_local')

# Row 2: Wrong
ax = axes[1, 0]
ax.scatter(sl_w[:, 0], sl_w[:, 1], s=100, c='indianred', alpha=0.7,
          edgecolors='darkred', linewidths=0.5)
ax.set_aspect('equal'); ax.grid(True, alpha=0.2)
ax.set_title(f'Raw Positions ({len(sl_w)} atoms)', fontsize=12)
ax.set_ylabel('WRONG dt=0.02\nY (σ)', fontsize=13, color='red', fontweight='bold')
ax.set_xlabel('X (σ)', fontsize=11)

cn_w = coordination_number(sl_w, r_cut=1.5)
ax = axes[1, 1]
sc3 = ax.scatter(sl_w[:, 0], sl_w[:, 1], s=100, c=cn_w, cmap='RdYlGn',
                vmin=1, vmax=8, alpha=0.8, edgecolors='k', linewidths=0.3)
ax.set_aspect('equal'); ax.grid(True, alpha=0.2)
ax.set_title(f'Coordination Number\nMean={cn_w.mean():.2f}, Std={cn_w.std():.2f}', fontsize=12)
ax.set_xlabel('X (σ)', fontsize=11)
plt.colorbar(sc3, ax=ax, shrink=0.7, label='CN')

ld_w = local_density(sl_w, r_probe=1.8)
ax = axes[1, 2]
sc4 = ax.scatter(sl_w[:, 0], sl_w[:, 1], s=100, c=ld_w, cmap='viridis',
                vmin=0.3, vmax=1.5, alpha=0.8, edgecolors='k', linewidths=0.3)
ax.set_aspect('equal'); ax.grid(True, alpha=0.2)
ax.set_title(f'Local Density (r<1.8σ)\nMean={ld_w.mean():.3f}', fontsize=12)
ax.set_xlabel('X (σ)', fontsize=11)
plt.colorbar(sc4, ax=ax, shrink=0.7, label='ρ_local')

fig.suptitle('Final Frame Diagnostic: How to Tell Wrong from Correct in XY Slice\n'
             '(Dense LJ Fluid ρ*=0.8442, T*=0.75)',
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/dense_diagnostic_panel.png", dpi=150)
plt.close()
print("  -> dense_diagnostic_panel.png")

# ============================================================
# Plot 3: Quantitative histograms side by side
# ============================================================
print("Generating quantitative histogram comparison...")

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# NN distance
nn_c = nearest_neighbor_dist(sl_c)
nn_w = nearest_neighbor_dist(sl_w)
axes[0].hist(nn_c, bins=40, alpha=0.6, color='steelblue', density=True, label=f'Correct (mean={nn_c.mean():.3f})')
axes[0].hist(nn_w, bins=40, alpha=0.6, color='indianred', density=True, label=f'Wrong (mean={nn_w.mean():.3f})')
axes[0].axvline(x=2**(1/6), color='green', linestyle='--', linewidth=2, label=f'LJ min={2**(1/6):.3f}σ')
axes[0].set_xlabel('Nearest-Neighbor Distance (σ)', fontsize=12)
axes[0].set_ylabel('Probability Density', fontsize=12)
axes[0].set_title('NN Distance Distribution', fontsize=13)
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3)

# Coordination number
axes[1].hist(cn_c, bins=range(0, 12), alpha=0.6, color='steelblue', density=True,
            label=f'Correct (mean={cn_c.mean():.2f})')
axes[1].hist(cn_w, bins=range(0, 12), alpha=0.6, color='indianred', density=True,
            label=f'Wrong (mean={cn_w.mean():.2f})')
axes[1].set_xlabel('Coordination Number', fontsize=12)
axes[1].set_ylabel('Probability', fontsize=12)
axes[1].set_title('Coordination Number Distribution', fontsize=13)
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)

# Local density
axes[2].hist(ld_c, bins=30, alpha=0.6, color='steelblue', density=True,
            label=f'Correct (mean={ld_c.mean():.3f})')
axes[2].hist(ld_w, bins=30, alpha=0.6, color='indianred', density=True,
            label=f'Wrong (mean={ld_w.mean():.3f})')
axes[2].set_xlabel('Local Density', fontsize=12)
axes[2].set_ylabel('Probability Density', fontsize=12)
axes[2].set_title('Local Density Distribution', fontsize=13)
axes[2].legend(fontsize=9)
axes[2].grid(True, alpha=0.3)

fig.suptitle('Quantitative Structural Metrics: Correct vs Wrong Timestep\n'
             'Wrong timestep → lower coordination, broader distributions',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/dense_histograms.png", dpi=150)
plt.close()
print("  -> dense_histograms.png")

# ============================================================
# Plot 4: Draw "bonds" (lines between close neighbors) to show structure
# ============================================================
print("Generating bond-network visualization...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))

r_bond = 1.3  # draw bonds shorter than this

for ax, sl, label, color, ecolor in [
    (ax1, sl_c, 'CORRECT dt=0.005', 'steelblue', 'navy'),
    (ax2, sl_w, 'WRONG dt=0.02', 'indianred', 'darkred')
]:
    # Draw bonds
    d = cdist(sl[:, :2], sl[:, :2])
    bond_count = 0
    for i in range(len(sl)):
        for j in range(i+1, len(sl)):
            if d[i, j] < r_bond:
                ax.plot([sl[i, 0], sl[j, 0]], [sl[i, 1], sl[j, 1]],
                       color='gray', linewidth=0.4, alpha=0.4)
                bond_count += 1

    # Draw atoms on top
    ax.scatter(sl[:, 0], sl[:, 1], s=80, c=color, alpha=0.8,
              edgecolors=ecolor, linewidths=0.5, zorder=5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2)
    ax.set_xlabel('X (σ)', fontsize=12)
    ax.set_ylabel('Y (σ)', fontsize=12)
    ax.set_title(f'{label}\n{len(sl)} atoms, {bond_count} bonds (r < {r_bond}σ)',
                fontsize=13, fontweight='bold')

fig.suptitle('Bond Network in XY Slice (Dense LJ Fluid ρ*=0.8442)\n'
             'Lines connect atoms closer than 1.3σ — denser network = more liquid-like ordering',
             fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/dense_bond_network.png", dpi=150)
plt.close()
print("  -> dense_bond_network.png")

print("\nAll dense comparison plots saved!")
