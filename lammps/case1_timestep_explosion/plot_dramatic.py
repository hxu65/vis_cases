#!/usr/bin/env python3
"""Dramatic XY slice comparison showing progressive structural breakdown."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist

BASE = "/mnt/common/hxu40/File/lammps_cases/case1_timestep_explosion"

def read_all_frames(filename):
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
        for j in range(i + 2, min(i + 2 + natoms, len(lines))):
            parts = lines[j].split()
            if len(parts) >= 4:
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        if coords:
            frames.append(np.array(coords))
        i += natoms + 2
    return frames

def xy_slice(coords, z_frac=0.1):
    zmid = np.median(coords[:, 2])
    zrange = coords[:, 2].max() - coords[:, 2].min()
    if zrange == 0:
        return coords
    z_lo = zmid - z_frac * zrange
    z_hi = zmid + z_frac * zrange
    mask = (coords[:, 2] >= z_lo) & (coords[:, 2] <= z_hi)
    return coords[mask]

def coordination_number(coords, r_cut=1.5):
    if len(coords) < 2:
        return np.zeros(len(coords))
    d = cdist(coords[:, :2], coords[:, :2])
    np.fill_diagonal(d, np.inf)
    return np.sum(d < r_cut, axis=1)

def draw_bonds(ax, sl, r_bond=1.3):
    d = cdist(sl[:, :2], sl[:, :2])
    count = 0
    for i in range(len(sl)):
        for j in range(i+1, len(sl)):
            if d[i, j] < r_bond:
                ax.plot([sl[i, 0], sl[j, 0]], [sl[i, 1], sl[j, 1]],
                       color='gray', linewidth=0.5, alpha=0.35)
                count += 1
    return count

# Load all trajectories
print("Loading trajectories...")
frames_correct = read_all_frames(f"{BASE}/trajectory_dense_correct.xyz")
frames_moderate = read_all_frames(f"{BASE}/trajectory_dense_moderate_wrong.xyz")
frames_bad = read_all_frames(f"{BASE}/trajectory_dense_bad.xyz")

print(f"  Correct (dt=0.005): {len(frames_correct)} frames")
print(f"  Moderate (dt=0.025): {len(frames_moderate)} frames")
print(f"  Bad (dt=0.03): {len(frames_bad)} frames")

# ============================================================
# PLOT 1: Progressive degradation — 4 regimes side by side
#         Bond network + atom coloring by coordination number
# ============================================================
print("\nGenerating progressive degradation plot...")

# Pick final frames from each
datasets = [
    (frames_correct[-1], 'CORRECT\ndt=0.005', 'steelblue', 'darkblue'),
    (frames_correct[0], 'CORRECT\n(start)', 'mediumseagreen', 'darkgreen'),
    (frames_moderate[-1], 'MODERATE WRONG\ndt=0.025 (step 1550)', 'orange', 'darkorange'),
    (frames_bad[-1], 'VERY WRONG\ndt=0.03 (step 110, pre-crash)', 'red', 'darkred'),
]

fig, axes = plt.subplots(1, 4, figsize=(28, 7))

for idx, (coords, label, color, ecolor) in enumerate(datasets):
    ax = axes[idx]
    sl = xy_slice(coords, z_frac=0.1)

    if len(sl) < 2:
        ax.text(0.5, 0.5, 'Not enough atoms\nin slice', transform=ax.transAxes,
               ha='center', fontsize=14)
        ax.set_title(label, fontsize=13, fontweight='bold')
        continue

    cn = coordination_number(sl, r_cut=1.5)
    bonds = draw_bonds(ax, sl, r_bond=1.3)

    sc = ax.scatter(sl[:, 0], sl[:, 1], s=120, c=cn, cmap='RdYlGn',
                   vmin=0, vmax=9, alpha=0.85, edgecolors='k', linewidths=0.4, zorder=5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.15)
    ax.set_xlabel('X (σ)', fontsize=12)
    if idx == 0:
        ax.set_ylabel('Y (σ)', fontsize=12)
    ax.set_title(f'{label}\n{len(sl)} atoms, {bonds} bonds\nMean CN={cn.mean():.2f}',
                fontsize=12, fontweight='bold')

plt.colorbar(sc, ax=axes, shrink=0.7, label='Coordination Number (r < 1.5σ)', pad=0.02)

fig.suptitle('Progressive Structural Breakdown from Wrong Timestep\n'
             '(Dense LJ Fluid ρ*=0.8442 — XY Slice with Bond Network)\n'
             'Green atoms = well-coordinated (liquid-like)  |  Red atoms = under-coordinated (gas-like)',
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/progressive_breakdown.png", dpi=150, bbox_inches='tight')
plt.close()
print("  -> progressive_breakdown.png")

# ============================================================
# PLOT 2: Time series for dt=0.025 — watch structure degrade
# ============================================================
print("Generating time-lapse degradation (dt=0.025)...")

n = len(frames_moderate)
pick = [0, n//6, n//3, n//2, 2*n//3, n-1]

fig, axes = plt.subplots(2, len(pick), figsize=(6*len(pick), 12))

for col, fidx in enumerate(pick):
    coords = frames_moderate[fidx]
    sl = xy_slice(coords, z_frac=0.1)

    if len(sl) < 2:
        for row in range(2):
            axes[row, col].text(0.5, 0.5, 'N/A', transform=axes[row, col].transAxes, ha='center')
        continue

    cn = coordination_number(sl, r_cut=1.5)

    # Top: bond network
    ax = axes[0, col]
    bonds = draw_bonds(ax, sl, r_bond=1.3)
    sc = ax.scatter(sl[:, 0], sl[:, 1], s=90, c=cn, cmap='RdYlGn',
                   vmin=0, vmax=9, alpha=0.85, edgecolors='k', linewidths=0.3, zorder=5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.15)
    step = fidx * 50  # dump every 50 steps
    ax.set_title(f'Step {step}\n{len(sl)} atoms, {bonds} bonds\nMean CN={cn.mean():.2f}',
                fontsize=11, fontweight='bold')
    if col == 0:
        ax.set_ylabel('Bond Network\nY (σ)', fontsize=12)

    # Bottom: atom size proportional to NN distance (big = isolated)
    ax = axes[1, col]
    nn = np.min(cdist(sl[:, :2], sl[:, :2]) + np.eye(len(sl))*999, axis=1)
    sizes = (nn / nn.min()) ** 3 * 40  # exaggerate size difference
    sc2 = ax.scatter(sl[:, 0], sl[:, 1], s=sizes, c=nn, cmap='hot_r',
                    vmin=0.85, vmax=1.8, alpha=0.8, edgecolors='k', linewidths=0.3)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.15)
    ax.set_xlabel('X (σ)', fontsize=11)
    if col == 0:
        ax.set_ylabel('NN Distance\nY (σ)', fontsize=12)
    ax.set_title(f'Mean NN={nn.mean():.3f}σ', fontsize=10)

fig.suptitle('Time-Lapse: Structure Degradation with Wrong Timestep (dt=0.025)\n'
             'Top: Bond network (green=coordinated, red=isolated)  |  '
             'Bottom: Atom size ∝ isolation (big dark = far from neighbors)',
             fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/timelapse_degradation.png", dpi=150, bbox_inches='tight')
plt.close()
print("  -> timelapse_degradation.png")

# ============================================================
# PLOT 3: dt=0.03 — frame-by-frame explosion sequence
# ============================================================
print("Generating explosion sequence (dt=0.03)...")

n_bad = len(frames_bad)
fig, axes = plt.subplots(2, min(n_bad, 6), figsize=(6*min(n_bad,6), 12))
if min(n_bad, 6) == 1:
    axes = axes.reshape(2, 1)

pick_bad = np.linspace(0, n_bad-1, min(n_bad, 6), dtype=int)

for col, fidx in enumerate(pick_bad):
    coords = frames_bad[fidx]
    sl = xy_slice(coords, z_frac=0.12)

    if len(sl) < 2:
        for row in range(2):
            axes[row, col].text(0.5, 0.5, 'N/A', transform=axes[row, col].transAxes, ha='center')
        continue

    cn = coordination_number(sl, r_cut=1.5)

    # Top: bond network colored by CN
    ax = axes[0, col]
    bonds = draw_bonds(ax, sl, r_bond=1.3)
    sc = ax.scatter(sl[:, 0], sl[:, 1], s=100, c=cn, cmap='RdYlGn',
                   vmin=0, vmax=9, alpha=0.85, edgecolors='k', linewidths=0.3, zorder=5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.15)
    step = fidx * 10
    ax.set_title(f'Step {step}\n{bonds} bonds, CN={cn.mean():.2f}', fontsize=11, fontweight='bold')
    if col == 0:
        ax.set_ylabel('Bond Network\nY (σ)', fontsize=12)

    # Bottom: raw positions with displacement vectors from initial
    ax = axes[1, col]
    sl0 = xy_slice(frames_bad[0], z_frac=0.12)
    # Just show positions colored by how many atoms are in the slice
    nn = np.min(cdist(sl[:, :2], sl[:, :2]) + np.eye(len(sl))*999, axis=1)
    sc2 = ax.scatter(sl[:, 0], sl[:, 1], s=100, c=nn, cmap='hot_r',
                    vmin=0.85, vmax=2.0, alpha=0.8, edgecolors='k', linewidths=0.3)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.15)
    ax.set_xlabel('X (σ)', fontsize=11)
    if col == 0:
        ax.set_ylabel('NN Distance\nY (σ)', fontsize=12)

fig.suptitle('Explosion Sequence (dt=0.03): Frame-by-Frame Structural Collapse\n'
             'Watch bonds disappear and atoms become isolated before crash at step ~110',
             fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/explosion_sequence_bonds.png", dpi=150, bbox_inches='tight')
plt.close()
print("  -> explosion_sequence_bonds.png")

# ============================================================
# PLOT 4: Summary statistics over time
# ============================================================
print("Generating statistics over time...")

def compute_stats_over_time(frames, z_frac=0.1, r_cut=1.5, r_bond=1.3):
    mean_cn, mean_nn, bond_counts, n_atoms = [], [], [], []
    for coords in frames:
        sl = xy_slice(coords, z_frac=z_frac)
        if len(sl) < 2:
            mean_cn.append(0); mean_nn.append(0); bond_counts.append(0); n_atoms.append(len(sl))
            continue
        cn = coordination_number(sl, r_cut=r_cut)
        nn = np.min(cdist(sl[:, :2], sl[:, :2]) + np.eye(len(sl))*999, axis=1)
        d = cdist(sl[:, :2], sl[:, :2])
        np.fill_diagonal(d, np.inf)
        bonds = np.sum(d < r_bond) // 2
        mean_cn.append(cn.mean())
        mean_nn.append(nn.mean())
        bond_counts.append(bonds)
        n_atoms.append(len(sl))
    return mean_cn, mean_nn, bond_counts, n_atoms

print("  Computing correct stats...")
cn_c, nn_c, bc_c, na_c = compute_stats_over_time(frames_correct)
print("  Computing moderate wrong stats...")
cn_m, nn_m, bc_m, na_m = compute_stats_over_time(frames_moderate)
print("  Computing bad stats...")
cn_b, nn_b, bc_b, na_b = compute_stats_over_time(frames_bad)

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Normalize x-axis to fraction of trajectory
x_c = np.linspace(0, 1, len(cn_c))
x_m = np.linspace(0, 1, len(cn_m))
x_b = np.linspace(0, 1, len(cn_b))

axes[0,0].plot(x_c, cn_c, 'b-', linewidth=1.5, label='Correct dt=0.005')
axes[0,0].plot(x_m, cn_m, 'orange', linewidth=1.5, label='Wrong dt=0.025')
axes[0,0].plot(x_b, cn_b, 'r-', linewidth=2, label='Bad dt=0.03')
axes[0,0].set_ylabel('Mean Coordination Number', fontsize=12)
axes[0,0].set_title('Coordination Number Over Time', fontsize=13)
axes[0,0].legend(fontsize=10)
axes[0,0].grid(True, alpha=0.3)

axes[0,1].plot(x_c, nn_c, 'b-', linewidth=1.5, label='Correct dt=0.005')
axes[0,1].plot(x_m, nn_m, 'orange', linewidth=1.5, label='Wrong dt=0.025')
axes[0,1].plot(x_b, nn_b, 'r-', linewidth=2, label='Bad dt=0.03')
axes[0,1].set_ylabel('Mean NN Distance (σ)', fontsize=12)
axes[0,1].set_title('Nearest-Neighbor Distance Over Time', fontsize=13)
axes[0,1].legend(fontsize=10)
axes[0,1].grid(True, alpha=0.3)

axes[1,0].plot(x_c, bc_c, 'b-', linewidth=1.5, label='Correct dt=0.005')
axes[1,0].plot(x_m, bc_m, 'orange', linewidth=1.5, label='Wrong dt=0.025')
axes[1,0].plot(x_b, bc_b, 'r-', linewidth=2, label='Bad dt=0.03')
axes[1,0].set_xlabel('Fraction of Trajectory', fontsize=12)
axes[1,0].set_ylabel('Bond Count (r < 1.3σ)', fontsize=12)
axes[1,0].set_title('Number of Bonds Over Time', fontsize=13)
axes[1,0].legend(fontsize=10)
axes[1,0].grid(True, alpha=0.3)

axes[1,1].plot(x_c, na_c, 'b-', linewidth=1.5, label='Correct dt=0.005')
axes[1,1].plot(x_m, na_m, 'orange', linewidth=1.5, label='Wrong dt=0.025')
axes[1,1].plot(x_b, na_b, 'r-', linewidth=2, label='Bad dt=0.03')
axes[1,1].set_xlabel('Fraction of Trajectory', fontsize=12)
axes[1,1].set_ylabel('Atoms in Slice', fontsize=12)
axes[1,1].set_title('Atoms in XY Slice Over Time', fontsize=13)
axes[1,1].legend(fontsize=10)
axes[1,1].grid(True, alpha=0.3)

fig.suptitle('Structural Metrics Over Time: Correct vs Wrong Timestep\n'
             'Wrong timestep → declining coordination, fewer bonds, increasing NN distance',
             fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/structural_metrics_time.png", dpi=150)
plt.close()
print("  -> structural_metrics_time.png")

print("\nAll dramatic comparison plots saved!")
