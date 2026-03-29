#!/usr/bin/env python3
"""XY slice comparison: stable vs explosion at multiple timeframes."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE = "/mnt/common/hxu40/File/lammps_cases/case1_timestep_explosion"

def read_all_frames(filename):
    """Read all frames from XYZ file, return list of coordinate arrays."""
    frames = []
    with open(filename) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        try:
            natoms = int(lines[i].strip())
        except:
            break
        coords = []
        for j in range(i + 2, i + 2 + natoms):
            parts = lines[j].split()
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        frames.append(np.array(coords))
        i += natoms + 2
    return frames

def xy_slice(coords, z_frac=0.1):
    """Extract a thin XY slab around the Z midplane."""
    zmid = np.median(coords[:, 2])
    zrange = coords[:, 2].max() - coords[:, 2].min()
    if zrange == 0:
        return coords
    z_lo = zmid - z_frac * zrange
    z_hi = zmid + z_frac * zrange
    mask = (coords[:, 2] >= z_lo) & (coords[:, 2] <= z_hi)
    return coords[mask]

def compute_pair_distances(coords, rmax=3.0):
    """Compute all pairwise distances < rmax for a 2D slice."""
    from scipy.spatial.distance import pdist
    dists = pdist(coords[:, :2])  # 2D distances in XY
    return dists[dists < rmax]

# ============================================================
# Load trajectories
# ============================================================
print("Loading trajectories...")
frames_correct = read_all_frames(f"{BASE}/trajectory_correct_v2.xyz")
frames_wrong = read_all_frames(f"{BASE}/trajectory_wrong_v2.xyz")
frames_explode = read_all_frames(f"{BASE}/trajectory.xyz")

print(f"  Correct: {len(frames_correct)} frames")
print(f"  Wrong (dt=0.02): {len(frames_wrong)} frames")
print(f"  Explosion (dt=0.05): {len(frames_explode)} frames")

# ============================================================
# Plot 1: Multi-frame XY slice evolution - Correct vs Wrong
# ============================================================
print("Generating multi-frame XY slice comparison...")

# Pick frames at 0%, 25%, 50%, 75%, 100% of trajectory
pcts = [0, 0.25, 0.5, 0.75, 1.0]
fig, axes = plt.subplots(2, 5, figsize=(25, 10))

for col, pct in enumerate(pcts):
    idx_c = min(int(pct * (len(frames_correct) - 1)), len(frames_correct) - 1)
    idx_w = min(int(pct * (len(frames_wrong) - 1)), len(frames_wrong) - 1)

    slice_c = xy_slice(frames_correct[idx_c])
    slice_w = xy_slice(frames_wrong[idx_w])

    # Correct (top row)
    ax = axes[0, col]
    ax.scatter(slice_c[:, 0], slice_c[:, 1], s=60, c='steelblue',
               alpha=0.7, edgecolors='darkblue', linewidths=0.3)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2)
    ax.set_title(f'Frame {idx_c} ({int(pct*100)}%)', fontsize=11)
    if col == 0:
        ax.set_ylabel('Correct dt=0.005\nY', fontsize=12, color='blue')

    # Wrong (bottom row)
    ax = axes[1, col]
    ax.scatter(slice_w[:, 0], slice_w[:, 1], s=60, c='indianred',
               alpha=0.7, edgecolors='darkred', linewidths=0.3)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2)
    ax.set_xlabel('X', fontsize=11)
    if col == 0:
        ax.set_ylabel('Wrong dt=0.02\nY', fontsize=12, color='red')

fig.suptitle('XY Slice Evolution: Correct vs Wrong Timestep\n(Thin slab at Z midplane, after equilibration)',
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/xy_slice_evolution.png", dpi=150)
plt.close()
print("  -> xy_slice_evolution.png")

# ============================================================
# Plot 2: Explosion sequence - frame by frame
# ============================================================
print("Generating explosion XY slice sequence...")

n_exp = len(frames_explode)
fig, axes = plt.subplots(1, n_exp, figsize=(5 * n_exp, 5))
if n_exp == 1:
    axes = [axes]

for i, frame in enumerate(frames_explode):
    ax = axes[i]
    sl = xy_slice(frame, z_frac=0.15)
    # Color by distance from center
    cx, cy = np.mean(frame[:, 0]), np.mean(frame[:, 1])
    if len(sl) > 0:
        dist = np.sqrt((sl[:, 0] - cx)**2 + (sl[:, 1] - cy)**2)
        sc = ax.scatter(sl[:, 0], sl[:, 1], s=50, c=dist,
                       cmap='hot_r', alpha=0.8, edgecolors='k', linewidths=0.2)
    ax.set_aspect('equal')
    ax.set_title(f'Step {i * 10}\n({len(sl)} atoms in slice)', fontsize=11)
    ax.grid(True, alpha=0.2)
    ax.set_xlabel('X', fontsize=10)
    if i == 0:
        ax.set_ylabel('Y', fontsize=10)

fig.suptitle('Explosion Sequence (dt=0.05): XY Slices Every 10 Steps\n(Color = distance from center)',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/xy_slice_explosion_sequence.png", dpi=150)
plt.close()
print("  -> xy_slice_explosion_sequence.png")

# ============================================================
# Plot 3: Nearest-neighbor distance distribution
# ============================================================
print("Generating nearest-neighbor analysis...")
from scipy.spatial.distance import cdist

def nearest_neighbor_dists(coords):
    """Compute nearest-neighbor distance for each atom."""
    if len(coords) < 2:
        return np.array([])
    d = cdist(coords, coords)
    np.fill_diagonal(d, np.inf)
    return np.min(d, axis=1)

# Compare last frames
nn_correct = nearest_neighbor_dists(frames_correct[-1])
nn_wrong = nearest_neighbor_dists(frames_wrong[-1])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

bins = np.linspace(0.5, 2.5, 80)
ax1.hist(nn_correct, bins=bins, color='steelblue', alpha=0.7, density=True, label='Correct dt=0.005')
ax1.hist(nn_wrong, bins=bins, color='indianred', alpha=0.7, density=True, label='Wrong dt=0.02')
ax1.axvline(x=2**(1/6), color='green', linestyle='--', linewidth=2,
            label=f'LJ min (r=2^(1/6)σ ≈ {2**(1/6):.3f})')
ax1.set_xlabel('Nearest-Neighbor Distance (σ)', fontsize=12)
ax1.set_ylabel('Probability Density', fontsize=12)
ax1.set_title('Nearest-Neighbor Distance Distribution', fontsize=13)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Mean squared displacement from initial positions
msd_correct, msd_wrong = [], []
ref_c = frames_correct[0]
ref_w = frames_wrong[0]
for i in range(len(frames_correct)):
    diff = frames_correct[i] - ref_c
    msd_correct.append(np.mean(np.sum(diff**2, axis=1)))
for i in range(len(frames_wrong)):
    diff = frames_wrong[i] - ref_w
    msd_wrong.append(np.mean(np.sum(diff**2, axis=1)))

ax2.plot(range(len(msd_correct)), msd_correct, 'b-', linewidth=1.5, label='Correct dt=0.005')
ax2.plot(range(len(msd_wrong)), msd_wrong, 'r-', linewidth=1.5, label='Wrong dt=0.02')
ax2.set_xlabel('Frame', fontsize=12)
ax2.set_ylabel('MSD (σ²)', fontsize=12)
ax2.set_title('Mean Squared Displacement from Initial Config', fontsize=13)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{BASE}/structural_analysis.png", dpi=150)
plt.close()
print("  -> structural_analysis.png")

# ============================================================
# Plot 4: Side-by-side final snapshot with Voronoi-like coloring
# ============================================================
print("Generating final snapshot comparison...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

for ax, frames, label, cmap_name in [
    (ax1, frames_correct, 'Correct (dt=0.005)', 'Blues'),
    (ax2, frames_wrong, 'Wrong (dt=0.02)', 'Reds')
]:
    coords = frames[-1]
    sl = xy_slice(coords, z_frac=0.12)
    nn = nearest_neighbor_dists(sl)

    sc = ax.scatter(sl[:, 0], sl[:, 1], s=100, c=nn, cmap=cmap_name,
                   alpha=0.8, edgecolors='k', linewidths=0.3,
                   vmin=0.9, vmax=1.8)
    plt.colorbar(sc, ax=ax, label='NN distance (σ)', shrink=0.8)
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_title(f'{label}\nMean NN dist = {nn.mean():.3f}σ', fontsize=13)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2)

fig.suptitle('Final Configuration: Atoms colored by Nearest-Neighbor Distance\n'
             '(More uniform color = more ordered structure)',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/xy_slice_nn_colored.png", dpi=150)
plt.close()
print("  -> xy_slice_nn_colored.png")

print("\nAll XY slice plots saved!")
