#!/usr/bin/env python3
"""Compare correct vs wrong timestep for Case 1."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE = "/mnt/common/hxu40/File/lammps_cases/case1_timestep_explosion"

def parse_thermo(filename):
    steps, temps, etotals = [], [], []
    with open(filename) as f:
        recording = False
        for line in f:
            if line.strip().startswith("Step"):
                recording = True
                continue
            if recording:
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        steps.append(float(parts[0]))
                        temps.append(float(parts[1]))
                        etotals.append(float(parts[4]))
                    except ValueError:
                        recording = False
    return np.array(steps), np.array(temps), np.array(etotals)

# Parse all three logs
s_correct, t_correct, e_correct = parse_thermo(f"{BASE}/log_correct.out")
s_mild, t_mild, e_mild = parse_thermo(f"{BASE}/log_mild.out")
s_explode, t_explode, e_explode = parse_thermo(f"{BASE}/log.out")

# ============================================================
# Plot 1: Correct vs Mild (dt=0.005 vs dt=0.01)
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(16, 10))

# Temperature comparison
axes[0, 0].plot(s_correct, t_correct, 'b-', linewidth=0.8, alpha=0.7, label='dt=0.005 (correct)')
axes[0, 0].set_ylabel('Temperature (T*)', fontsize=12)
axes[0, 0].set_title('Correct Timestep (dt=0.005)', fontsize=13, color='blue')
axes[0, 0].axhline(y=0.75, color='gray', linestyle='--', alpha=0.5, label='Initial T*')
axes[0, 0].legend(fontsize=10)
axes[0, 0].grid(True, alpha=0.3)
axes[0, 0].set_ylim(0.4, 1.2)

axes[0, 1].plot(s_mild, t_mild, 'r-', linewidth=0.8, alpha=0.7, label='dt=0.01 (2x too large)')
axes[0, 1].set_ylabel('Temperature (T*)', fontsize=12)
axes[0, 1].set_title('Wrong Timestep (dt=0.01)', fontsize=13, color='red')
axes[0, 1].axhline(y=0.75, color='gray', linestyle='--', alpha=0.5, label='Initial T*')
axes[0, 1].legend(fontsize=10)
axes[0, 1].grid(True, alpha=0.3)
axes[0, 1].set_ylim(0.4, 1.2)

# Energy conservation comparison
e_drift_correct = e_correct - e_correct[0]
e_drift_mild = e_mild - e_mild[0]

axes[1, 0].plot(s_correct, e_drift_correct, 'b-', linewidth=0.8)
axes[1, 0].set_xlabel('Timestep', fontsize=12)
axes[1, 0].set_ylabel('Energy Drift (E - E₀)', fontsize=12)
axes[1, 0].axhline(y=0, color='gray', linestyle='--', alpha=0.5)
axes[1, 0].set_title(f'Energy Drift: max |ΔE/E₀| = {np.max(np.abs(e_drift_correct))/abs(e_correct[0]):.2e}', fontsize=11)
axes[1, 0].grid(True, alpha=0.3)

axes[1, 1].plot(s_mild, e_drift_mild, 'r-', linewidth=0.8)
axes[1, 1].set_xlabel('Timestep', fontsize=12)
axes[1, 1].set_ylabel('Energy Drift (E - E₀)', fontsize=12)
axes[1, 1].axhline(y=0, color='gray', linestyle='--', alpha=0.5)
axes[1, 1].set_title(f'Energy Drift: max |ΔE/E₀| = {np.max(np.abs(e_drift_mild))/abs(e_mild[0]):.2e}', fontsize=11)
axes[1, 1].grid(True, alpha=0.3)

fig.suptitle('Case 1: Correct vs Wrong Timestep — NVE Energy Conservation', fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/correct_vs_wrong_comparison.png", dpi=150)
plt.close()
print("-> correct_vs_wrong_comparison.png")

# ============================================================
# Plot 2: All three on same axes for direct comparison
# ============================================================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), sharex=False)

# Temperature overlay (exclude explosion since it's on a totally different scale)
ax1.plot(s_correct, t_correct, 'b-', linewidth=0.8, alpha=0.7, label='dt=0.005 (correct)')
ax1.plot(s_mild, t_mild, 'r-', linewidth=0.8, alpha=0.7, label='dt=0.01 (wrong)')
ax1.axhline(y=0.75, color='gray', linestyle='--', alpha=0.5, label='Initial T*=0.75')
ax1.set_ylabel('Temperature (T*)', fontsize=13)
ax1.set_title('Temperature Evolution: Correct vs Wrong Timestep', fontsize=14)
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3)

# Energy drift overlay
ax2.plot(s_correct, (e_correct - e_correct[0]) / abs(e_correct[0]) * 100,
         'b-', linewidth=1, label=f'dt=0.005: max drift = {np.max(np.abs(e_drift_correct))/abs(e_correct[0])*100:.2f}%')
ax2.plot(s_mild, (e_mild - e_mild[0]) / abs(e_mild[0]) * 100,
         'r-', linewidth=1, label=f'dt=0.01: max drift = {np.max(np.abs(e_drift_mild))/abs(e_mild[0])*100:.2f}%')
ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax2.set_xlabel('Timestep', fontsize=13)
ax2.set_ylabel('Energy Drift (%)', fontsize=13)
ax2.set_title('Total Energy Drift (NVE should conserve energy)', fontsize=14)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{BASE}/overlay_comparison.png", dpi=150)
plt.close()
print("-> overlay_comparison.png")

# ============================================================
# Plot 3: The catastrophic explosion on its own
# ============================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

ax1.plot(s_explode, t_explode, 'r-o', linewidth=2, markersize=8)
ax1.set_yscale('log')
ax1.set_xlabel('Timestep', fontsize=12)
ax1.set_ylabel('Temperature (T*)', fontsize=12)
ax1.set_title('dt=0.05 (10x wrong): Catastrophic Blowup', fontsize=13)
ax1.grid(True, alpha=0.3)
ax1.annotate(f'T = {t_explode[-1]:.0e}!', xy=(s_explode[-1], t_explode[-1]),
             fontsize=12, color='red', fontweight='bold',
             xytext=(-60, -20), textcoords='offset points')

# Correct for context
ax2.plot(s_correct[:20], t_correct[:20], 'b-o', linewidth=2, markersize=5)
ax2.set_xlabel('Timestep', fontsize=12)
ax2.set_ylabel('Temperature (T*)', fontsize=12)
ax2.set_title('dt=0.005 (correct): Stable', fontsize=13)
ax2.grid(True, alpha=0.3)
ax2.set_ylim(0, 1.5)

fig.suptitle('Case 1: What Explosion vs Stability Looks Like', fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/explosion_vs_stable.png", dpi=150)
plt.close()
print("-> explosion_vs_stable.png")

# ============================================================
# Plot 4: Atom snapshots - correct vs mild
# ============================================================
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
    # Return last frame for frame_idx == -1
    natoms = int(lines[last_start].strip())
    coords = []
    for j in range(last_start + 2, last_start + 2 + natoms):
        parts = lines[j].split()
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.array(coords)

coords_correct_last = read_xyz_frame(f"{BASE}/trajectory_correct.xyz", -1)
coords_mild_last = read_xyz_frame(f"{BASE}/trajectory_mild.xyz", -1)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
for ax, coords, label, color in [
    (ax1, coords_correct_last, 'Correct (dt=0.005)', 'blue'),
    (ax2, coords_mild_last, 'Wrong (dt=0.01)', 'red')
]:
    zmid = np.median(coords[:, 2])
    zrange = coords[:, 2].max() - coords[:, 2].min()
    mask = (coords[:, 2] >= zmid - 0.1*zrange) & (coords[:, 2] <= zmid + 0.1*zrange)
    sliced = coords[mask]
    ax.scatter(sliced[:, 0], sliced[:, 1], s=80, c=color,
               alpha=0.6, edgecolors='k', linewidths=0.3)
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_title(f'{label} ({len(sliced)} atoms in slice)', fontsize=13)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

fig.suptitle('Final Atom Positions: Correct vs Wrong Timestep (XY slice)', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/atoms_correct_vs_wrong.png", dpi=150)
plt.close()
print("-> atoms_correct_vs_wrong.png")

print("\nDone! All comparison plots saved.")
