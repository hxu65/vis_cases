#!/usr/bin/env python3
"""Proper comparison: equilibrated system, correct vs wrong timestep."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE = "/mnt/common/hxu40/File/lammps_cases/case1_timestep_explosion"

def parse_thermo_production(filename):
    """Parse only the SECOND run block (production NVE) from the log."""
    steps, temps, etotals = [], [], []
    run_count = 0
    recording = False
    with open(filename) as f:
        for line in f:
            if line.strip().startswith("Step"):
                run_count += 1
                recording = (run_count == 2)  # only 2nd run
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

s_c, t_c, e_c = parse_thermo_production(f"{BASE}/log_correct_v2.out")
s_w, t_w, e_w = parse_thermo_production(f"{BASE}/log_wrong_v2.out")

# ============================================================
# Main comparison figure
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(16, 10))

# Temperature
axes[0, 0].plot(s_c, t_c, 'b-', linewidth=0.8)
axes[0, 0].set_ylabel('Temperature (T*)', fontsize=12)
axes[0, 0].set_title('Correct: dt = 0.005', fontsize=14, color='blue')
axes[0, 0].grid(True, alpha=0.3)
axes[0, 0].set_ylim(0.85, 1.15)

axes[0, 1].plot(s_w, t_w, 'r-', linewidth=0.8)
axes[0, 1].set_ylabel('Temperature (T*)', fontsize=12)
axes[0, 1].set_title('Wrong: dt = 0.02 (4x too large)', fontsize=14, color='red')
axes[0, 1].grid(True, alpha=0.3)
axes[0, 1].set_ylim(0.85, 1.15)

# Energy drift
drift_c = (e_c - e_c[0]) / abs(e_c[0]) * 100
drift_w = (e_w - e_w[0]) / abs(e_w[0]) * 100

axes[1, 0].plot(s_c, drift_c, 'b-', linewidth=1)
axes[1, 0].axhline(y=0, color='gray', linestyle='--', alpha=0.5)
axes[1, 0].set_xlabel('Timestep', fontsize=12)
axes[1, 0].set_ylabel('Energy Drift (%)', fontsize=12)
axes[1, 0].set_title(f'Energy Conservation: max drift = {np.max(np.abs(drift_c)):.3f}%', fontsize=12)
axes[1, 0].grid(True, alpha=0.3)

axes[1, 1].plot(s_w, drift_w, 'r-', linewidth=1)
axes[1, 1].axhline(y=0, color='gray', linestyle='--', alpha=0.5)
axes[1, 1].set_xlabel('Timestep', fontsize=12)
axes[1, 1].set_ylabel('Energy Drift (%)', fontsize=12)
axes[1, 1].set_title(f'Energy Drift: max drift = {np.max(np.abs(drift_w)):.3f}%', fontsize=12)
axes[1, 1].grid(True, alpha=0.3)

fig.suptitle('Case 1: Good vs Bad Timestep (Pre-Equilibrated NVE Production)',
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/good_vs_bad_v2.png", dpi=150)
plt.close()
print("-> good_vs_bad_v2.png")

# ============================================================
# Overlay figure
# ============================================================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9))

ax1.plot(s_c, t_c, 'b-', linewidth=0.8, alpha=0.7, label='dt=0.005 (correct)')
ax1.plot(s_w, t_w, 'r-', linewidth=0.8, alpha=0.7, label='dt=0.02 (wrong)')
ax1.set_ylabel('Temperature (T*)', fontsize=13)
ax1.set_title('Temperature: Correct vs Wrong Timestep (after equilibration)', fontsize=14)
ax1.legend(fontsize=12)
ax1.grid(True, alpha=0.3)

ax2.plot(s_c, drift_c, 'b-', linewidth=1.2,
         label=f'dt=0.005: max drift = {np.max(np.abs(drift_c)):.3f}%')
ax2.plot(s_w, drift_w, 'r-', linewidth=1.2,
         label=f'dt=0.02: max drift = {np.max(np.abs(drift_w)):.3f}%')
ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax2.set_xlabel('Timestep', fontsize=13)
ax2.set_ylabel('Energy Drift (%)', fontsize=13)
ax2.set_title('NVE Energy Conservation (should be flat at 0%)', fontsize=14)
ax2.legend(fontsize=12)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{BASE}/overlay_v2.png", dpi=150)
plt.close()
print("-> overlay_v2.png")

# ============================================================
# Explosion vs stable (all 3 regimes)
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Parse explosion log
s_e, t_e = [], []
with open(f"{BASE}/log.out") as f:
    recording = False
    for line in f:
        if line.strip().startswith("Step"):
            recording = True
            continue
        if recording:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    s_e.append(float(parts[0]))
                    t_e.append(float(parts[1]))
                except ValueError:
                    recording = False

axes[0].plot(s_e, t_e, 'darkred', marker='o', linewidth=2, markersize=8)
axes[0].set_yscale('log')
axes[0].set_title('dt=0.05 (10x wrong)\nCATASTROPHIC', fontsize=13, color='darkred')
axes[0].set_xlabel('Step', fontsize=11)
axes[0].set_ylabel('Temperature', fontsize=11)
axes[0].grid(True, alpha=0.3)

axes[1].plot(s_w, t_w, 'r-', linewidth=1)
axes[1].set_title('dt=0.02 (4x wrong)\nSLOW DRIFT', fontsize=13, color='red')
axes[1].set_xlabel('Step', fontsize=11)
axes[1].set_ylabel('Temperature', fontsize=11)
axes[1].set_ylim(0.85, 1.15)
axes[1].grid(True, alpha=0.3)

axes[2].plot(s_c, t_c, 'b-', linewidth=1)
axes[2].set_title('dt=0.005 (correct)\nSTABLE', fontsize=13, color='blue')
axes[2].set_xlabel('Step', fontsize=11)
axes[2].set_ylabel('Temperature', fontsize=11)
axes[2].set_ylim(0.85, 1.15)
axes[2].grid(True, alpha=0.3)

fig.suptitle('Three Regimes of Timestep Error', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{BASE}/three_regimes.png", dpi=150)
plt.close()
print("-> three_regimes.png")

print("\nAll v2 plots saved!")
