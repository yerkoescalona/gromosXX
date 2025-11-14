"""
Example 4: Polymer Chain Dynamics
==================================

Demonstrates:
- Covalent bonds (harmonic potentials)
- Bond angle potentials
- Polymer flexibility
- Radius of gyration
- Entropic elasticity
- Thermal fluctuations

This simulation shows a flexible polymer chain with bonded interactions,
demonstrating how bonds and angles constrain molecular motion while
thermal energy drives conformational changes.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from core import Simulation2D, Renderer2D

# Simulation parameters
CHAIN_LENGTH = 20  # Number of monomers
BOX_SIZE = (6.0, 6.0)
TEMPERATURE = 300.0  # K
BOND_LENGTH = 0.15  # nm
BOND_K = 5000.0  # kJ/(mol·nm²) - stiff bond
ANGLE_K = 100.0  # kJ/(mol·rad²) - moderately stiff

print("=" * 70)
print("Educational 2D Simulation: Polymer Chain Dynamics")
print("=" * 70)
print(f"\nSimulation Parameters:")
print(f"  Chain length: {CHAIN_LENGTH} monomers")
print(f"  Box size: {BOX_SIZE[0]} x {BOX_SIZE[1]} nm")
print(f"  Temperature: {TEMPERATURE} K")
print(f"  Bond length: {BOND_LENGTH} nm")
print(f"  Bond stiffness: {BOND_K} kJ/(mol·nm²)")
print(f"  Angle stiffness: {ANGLE_K} kJ/(mol·rad²)")
print(f"  Periodic boundaries: Yes")
print()

# Create simulation
sim = Simulation2D(
    box_size=BOX_SIZE,
    periodic=True,
    integrator='velocity-verlet',
    dt=0.001  # Small timestep for stiff bonds
)

# Adjust nonbonded LJ parameters (softer for polymer)
sim.lj_params[0] = {
    'C6': 0.001,
    'C12': 0.000005,
    'sigma': 0.25,
    'epsilon': 0.3
}

print("Building polymer chain...")

# Create linear polymer chain starting from center
center = np.array(BOX_SIZE) / 2
current_pos = center.copy()
angle = 0.0

particle_indices = []

for i in range(CHAIN_LENGTH):
    # Position along chain with small random perturbation
    if i == 0:
        pos = current_pos
    else:
        # Random walk with preferred direction
        angle += np.random.uniform(-np.pi/4, np.pi/4)
        direction = np.array([np.cos(angle), np.sin(angle)])
        pos = current_pos + BOND_LENGTH * direction

    # Color coding: head (blue) -> middle (green) -> tail (red)
    if i < CHAIN_LENGTH / 3:
        color = 'blue'
    elif i < 2 * CHAIN_LENGTH / 3:
        color = 'green'
    else:
        color = 'red'

    idx = sim.add_particle(
        position=pos,
        mass=1.0,
        atom_type=0,
        color=color
    )
    particle_indices.append(idx)
    current_pos = pos

# Add bonds between consecutive monomers
print(f"Adding {CHAIN_LENGTH-1} bonds...")
for i in range(CHAIN_LENGTH - 1):
    sim.add_bond(
        i=particle_indices[i],
        j=particle_indices[i+1],
        k=BOND_K,
        r0=BOND_LENGTH
    )

# Add angle potentials for i-j-k triplets
print(f"Adding {CHAIN_LENGTH-2} angle potentials...")
for i in range(CHAIN_LENGTH - 2):
    sim.add_angle(
        i=particle_indices[i],
        j=particle_indices[i+1],
        k=particle_indices[i+2],
        k_angle=ANGLE_K,
        theta0=np.pi  # Prefer extended conformation
    )

print(f"\nPolymer chain complete:")
print(f"  {CHAIN_LENGTH} monomers")
print(f"  {len(sim.bonds)} bonds")
print(f"  {len(sim.angles)} angle potentials")

# Initialize velocities
print(f"\nInitializing velocities at {TEMPERATURE} K...")
sim.initialize_velocities(TEMPERATURE)
sim.set_temperature(TEMPERATURE, tau=0.1)

print(f"Initial temperature: {sim.temperature():.1f} K")

# Calculate initial radius of gyration
def calculate_radius_of_gyration(positions):
    """Calculate radius of gyration of polymer."""
    com = np.mean(positions, axis=0)
    distances_sq = np.sum((positions - com)**2, axis=1)
    return np.sqrt(np.mean(distances_sq))

def calculate_end_to_end_distance(positions):
    """Calculate end-to-end distance."""
    return np.linalg.norm(positions[-1] - positions[0])

# Create renderer with trails to show motion
renderer = Renderer2D(
    box_size=BOX_SIZE,
    figsize=(14, 6),
    show_forces=False,
    show_velocities=False,
    show_trails=True,
    trail_length=100
)

# Track conformational properties
rg_values = []
ree_values = []

# Phase 1: Equilibration
print("\n" + "=" * 70)
print("Phase 1: Equilibration (1000 steps)")
print("=" * 70)
print("\nObserve: Polymer explores conformational space")
print()

for step in range(1000):
    sim.integrate_step()

    # Calculate properties
    chain_positions = sim.get_positions()[:CHAIN_LENGTH]
    rg = calculate_radius_of_gyration(chain_positions)
    ree = calculate_end_to_end_distance(chain_positions)

    rg_values.append(rg)
    ree_values.append(ree)

    if step % 10 == 0:
        renderer.render(sim, show_bonds=True)

    if step % 200 == 0:
        ke, pe, total = sim.total_energy()
        print(f"Step {step:4d}: T = {sim.temperature():6.1f} K, "
              f"Rg = {rg:.3f} nm, Ree = {ree:.3f} nm, "
              f"E_bond = {sim.energies['bond'][-1]:6.2f} kJ/mol")

# Phase 2: Production run with analysis
print("\n" + "=" * 70)
print("Phase 2: Production Run (2000 steps)")
print("=" * 70)
print("\nCollecting conformational statistics...")
print()

for step in range(2000):
    sim.integrate_step()

    chain_positions = sim.get_positions()[:CHAIN_LENGTH]
    rg = calculate_radius_of_gyration(chain_positions)
    ree = calculate_end_to_end_distance(chain_positions)

    rg_values.append(rg)
    ree_values.append(ree)

    if step % 10 == 0:
        renderer.render(sim, show_bonds=True)

    if step % 400 == 0:
        ke, pe, total = sim.total_energy()
        print(f"Step {step:4d}: T = {sim.temperature():6.1f} K, "
              f"Rg = {rg:.3f} nm, Ree = {ree:.3f} nm")

print("\n" + "=" * 70)
print("Simulation Complete!")
print("=" * 70)

# Final analysis
chain_positions = sim.get_positions()[:CHAIN_LENGTH]
final_rg = calculate_radius_of_gyration(chain_positions)
final_ree = calculate_end_to_end_distance(chain_positions)

# Statistics
mean_rg = np.mean(rg_values)
std_rg = np.std(rg_values)
mean_ree = np.mean(ree_values)
std_ree = np.std(ree_values)

ke, pe, total = sim.total_energy()
bond_energy = sim.energies['bond'][-1]
angle_energy = sim.energies['angle'][-1]

print(f"\nFinal Statistics:")
print(f"  Temperature: {sim.temperature():.1f} K")
print(f"  Final Rg: {final_rg:.3f} nm")
print(f"  Final Ree: {final_ree:.3f} nm")
print(f"  Average Rg: {mean_rg:.3f} ± {std_rg:.3f} nm")
print(f"  Average Ree: {mean_ree:.3f} ± {std_ree:.3f} nm")
print(f"  Bond energy: {bond_energy:.2f} kJ/mol")
print(f"  Angle energy: {angle_energy:.2f} kJ/mol")

# Theoretical prediction for 2D random walk
# <Ree²> = N * b² for freely-jointed chain
# where N = number of bonds, b = bond length
expected_ree_sq = (CHAIN_LENGTH - 1) * BOND_LENGTH**2
expected_ree = np.sqrt(expected_ree_sq)

print(f"\nTheoretical Analysis (2D Random Walk):")
print(f"  Expected Ree (freely-jointed): {expected_ree:.3f} nm")
print(f"  Observed Ree: {mean_ree:.3f} nm")
print(f"  Ratio: {mean_ree/expected_ree:.2f} (>1 = extended, <1 = compact)")

# Save final frame
output_file = os.path.join(os.path.dirname(__file__), '../assets/polymer_chain.png')
renderer.save_frame(output_file)
print(f"\nSaved final frame to: {output_file}")

# Create energy summary
from core.renderer import plot_energy_summary

fig = plot_energy_summary(sim, figsize=(12, 8))
summary_file = os.path.join(os.path.dirname(__file__), '../assets/polymer_chain_energy.png')
fig.savefig(summary_file, dpi=150, bbox_inches='tight')
print(f"Saved energy summary to: {summary_file}")

# Create conformational analysis plot
fig_conf, axes = plt.subplots(2, 2, figsize=(12, 10))

times = np.arange(len(rg_values)) * sim.dt

# Radius of gyration
axes[0, 0].plot(times, rg_values, 'b-', linewidth=1, alpha=0.7)
axes[0, 0].axhline(mean_rg, color='r', linestyle='--', label=f'Mean: {mean_rg:.3f} nm')
axes[0, 0].set_xlabel('Time (ps)')
axes[0, 0].set_ylabel('Radius of Gyration (nm)')
axes[0, 0].set_title('Conformational Fluctuations: Rg')
axes[0, 0].legend()
axes[0, 0].grid(True, alpha=0.3)

# End-to-end distance
axes[0, 1].plot(times, ree_values, 'g-', linewidth=1, alpha=0.7)
axes[0, 1].axhline(mean_ree, color='r', linestyle='--', label=f'Mean: {mean_ree:.3f} nm')
axes[0, 1].axhline(expected_ree, color='orange', linestyle=':',
                   label=f'Expected (free): {expected_ree:.3f} nm')
axes[0, 1].set_xlabel('Time (ps)')
axes[0, 1].set_ylabel('End-to-End Distance (nm)')
axes[0, 1].set_title('End-to-End Distance')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.3)

# Histogram of Rg
axes[1, 0].hist(rg_values, bins=30, alpha=0.7, edgecolor='black')
axes[1, 0].axvline(mean_rg, color='r', linestyle='--', linewidth=2,
                   label=f'Mean: {mean_rg:.3f} nm')
axes[1, 0].set_xlabel('Radius of Gyration (nm)')
axes[1, 0].set_ylabel('Frequency')
axes[1, 0].set_title('Rg Distribution')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3)

# Histogram of Ree
axes[1, 1].hist(ree_values, bins=30, alpha=0.7, edgecolor='black', color='green')
axes[1, 1].axvline(mean_ree, color='r', linestyle='--', linewidth=2,
                   label=f'Mean: {mean_ree:.3f} nm')
axes[1, 1].axvline(expected_ree, color='orange', linestyle=':', linewidth=2,
                   label=f'Free chain: {expected_ree:.3f} nm')
axes[1, 1].set_xlabel('End-to-End Distance (nm)')
axes[1, 1].set_ylabel('Frequency')
axes[1, 1].set_title('Ree Distribution')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
conf_file = os.path.join(os.path.dirname(__file__), '../assets/polymer_chain_conformations.png')
fig_conf.savefig(conf_file, dpi=150, bbox_inches='tight')
print(f"Saved conformational analysis to: {conf_file}")

print("\nEducational Notes:")
print("  - Bonds constrain distances between neighbors")
print("  - Angles provide some stiffness to the chain")
print("  - Thermal fluctuations drive conformational changes")
print("  - Rg and Ree characterize polymer size/shape")
print("  - This is a 2D version of polymer physics!")

renderer.show()
