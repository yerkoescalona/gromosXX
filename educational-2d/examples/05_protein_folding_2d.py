"""
Example 5: Simple Protein Folding (2D Model)
=============================================

Demonstrates:
- Hydrophobic collapse
- Different atom types (hydrophobic vs hydrophilic)
- Secondary structure formation
- Folding funnel
- Native state stability

This simulation shows a simplified 2D protein model where hydrophobic
residues (represented by one color) prefer to cluster together,
driving the chain to fold into a compact structure.

This is inspired by the HP (Hydrophobic-Polar) protein model used in
computational biology education.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from core import Simulation2D, Renderer2D

# Simulation parameters
CHAIN_LENGTH = 25
BOX_SIZE = (8.0, 8.0)
TEMPERATURE = 200.0  # K (moderate for folding)
BOND_LENGTH = 0.15
BOND_K = 5000.0
ANGLE_K = 50.0  # More flexible than polymer

# Define a sequence: H = hydrophobic (1), P = polar (0)
# Design a sequence that should fold into a compact structure
# Using a pattern that creates a hydrophobic core
SEQUENCE = [
    1, 1, 0, 0, 1, 1, 1, 0, 0, 1,  # H H P P H H H P P H
    1, 1, 0, 0, 1, 1, 1, 0, 0, 1,  # H H P P H H H P P H
    1, 1, 0, 0, 1                   # H H P P H
]

print("=" * 70)
print("Educational 2D Simulation: Protein Folding")
print("=" * 70)
print(f"\nSimulation Parameters:")
print(f"  Chain length: {CHAIN_LENGTH} residues")
print(f"  Box size: {BOX_SIZE[0]} x {BOX_SIZE[1]} nm")
print(f"  Temperature: {TEMPERATURE} K")
print(f"  Model: Hydrophobic-Polar (HP) in 2D")
print()

# Create simulation
sim = Simulation2D(
    box_size=BOX_SIZE,
    periodic=False,  # No PBC for folding
    integrator='velocity-verlet',
    dt=0.001
)

# Setup different LJ parameters for different residue types
# Type 0 (P - Polar): Weak self-interaction
sim.lj_params[0] = {
    'C6': 0.0005,
    'C12': 0.000001,
    'sigma': 0.25,
    'epsilon': 0.2
}

# Type 1 (H - Hydrophobic): Strong self-interaction (drives collapse)
sim.lj_params[1] = {
    'C6': 0.004,     # Strong attraction between hydrophobic residues
    'C12': 0.00002,
    'sigma': 0.30,
    'epsilon': 1.0
}

# Type 2 (H-P interaction): Weak interaction (mixing rule will be applied)
sim.lj_params[2] = {
    'C6': 0.001,
    'C12': 0.000005,
    'sigma': 0.275,
    'epsilon': 0.4
}

print("Protein Sequence (H=Hydrophobic, P=Polar):")
seq_str = ''.join(['H' if s == 1 else 'P' for s in SEQUENCE])
print(f"  {seq_str}")
print(f"  Hydrophobic: {sum(SEQUENCE)} / {len(SEQUENCE)} residues")
print(f"  Polar: {len(SEQUENCE) - sum(SEQUENCE)} / {len(SEQUENCE)} residues")
print()

print("Building protein chain...")

# Start from center in extended conformation
center = np.array(BOX_SIZE) / 2
current_pos = center - np.array([CHAIN_LENGTH * BOND_LENGTH / 2, 0])

particle_indices = []

for i in range(CHAIN_LENGTH):
    # Initial extended conformation along x-axis
    pos = current_pos + np.array([i * BOND_LENGTH, 0])

    # Color based on residue type
    residue_type = SEQUENCE[i]
    if residue_type == 1:
        color = 'red'        # Hydrophobic (red)
        atom_type = 1
    else:
        color = 'lightblue'  # Polar (blue)
        atom_type = 0

    idx = sim.add_particle(
        position=pos,
        mass=1.0,
        atom_type=atom_type,
        color=color
    )
    particle_indices.append(idx)

# Add backbone bonds
print(f"Adding backbone bonds and angles...")
for i in range(CHAIN_LENGTH - 1):
    sim.add_bond(
        i=particle_indices[i],
        j=particle_indices[i+1],
        k=BOND_K,
        r0=BOND_LENGTH
    )

# Add backbone angles (more flexible for folding)
for i in range(CHAIN_LENGTH - 2):
    sim.add_angle(
        i=particle_indices[i],
        j=particle_indices[i+1],
        k=particle_indices[i+2],
        k_angle=ANGLE_K,
        theta0=np.pi  # Prefer extended, but flexible
    )

print(f"\nProtein structure:")
print(f"  {CHAIN_LENGTH} residues")
print(f"  {len(sim.bonds)} backbone bonds")
print(f"  {len(sim.angles)} backbone angles")

# Initialize velocities
print(f"\nInitializing velocities at {TEMPERATURE} K...")
sim.initialize_velocities(TEMPERATURE)
sim.set_temperature(TEMPERATURE, tau=0.2)

print(f"Initial temperature: {sim.temperature():.1f} K")

# Calculate properties
def calculate_radius_of_gyration(positions):
    com = np.mean(positions, axis=0)
    return np.sqrt(np.mean(np.sum((positions - com)**2, axis=1)))

def calculate_hydrophobic_contacts(positions, sequence, cutoff=0.4):
    """Count hydrophobic-hydrophobic contacts within cutoff distance."""
    contacts = 0
    n = len(positions)
    for i in range(n):
        if sequence[i] == 1:  # Hydrophobic
            for j in range(i+3, n):  # Skip nearest neighbors
                if sequence[j] == 1:  # Also hydrophobic
                    r = np.linalg.norm(positions[i] - positions[j])
                    if r < cutoff:
                        contacts += 1
    return contacts

# Create renderer with trails
renderer = Renderer2D(
    box_size=BOX_SIZE,
    figsize=(14, 6),
    show_forces=False,
    show_velocities=False,
    show_trails=True,
    trail_length=50
)

# Track folding progress
rg_values = []
contacts_values = []

# Phase 1: Initial collapse
print("\n" + "=" * 70)
print("Phase 1: Initial Folding (2000 steps)")
print("=" * 70)
print("\nObserve: Hydrophobic residues (red) begin to cluster")
print("         Polar residues (blue) stay on the outside")
print()

for step in range(2000):
    sim.integrate_step()

    chain_positions = sim.get_positions()[:CHAIN_LENGTH]
    rg = calculate_radius_of_gyration(chain_positions)
    contacts = calculate_hydrophobic_contacts(chain_positions, SEQUENCE)

    rg_values.append(rg)
    contacts_values.append(contacts)

    if step % 10 == 0:
        renderer.render(sim, show_bonds=True)

    if step % 400 == 0:
        ke, pe, total = sim.total_energy()
        print(f"Step {step:4d}: T = {sim.temperature():6.1f} K, "
              f"Rg = {rg:.3f} nm, H-contacts = {contacts:2d}, "
              f"E_pot = {pe:7.2f} kJ/mol")

# Phase 2: Continued equilibration
print("\n" + "=" * 70)
print("Phase 2: Equilibration (2000 steps)")
print("=" * 70)
print("\nWatching for stable folded state...")
print()

for step in range(2000):
    sim.integrate_step()

    chain_positions = sim.get_positions()[:CHAIN_LENGTH]
    rg = calculate_radius_of_gyration(chain_positions)
    contacts = calculate_hydrophobic_contacts(chain_positions, SEQUENCE)

    rg_values.append(rg)
    contacts_values.append(contacts)

    if step % 10 == 0:
        renderer.render(sim, show_bonds=True)

    if step % 400 == 0:
        ke, pe, total = sim.total_energy()
        print(f"Step {step:4d}: T = {sim.temperature():6.1f} K, "
              f"Rg = {rg:.3f} nm, H-contacts = {contacts:2d}")

print("\n" + "=" * 70)
print("Simulation Complete!")
print("=" * 70)

# Final analysis
chain_positions = sim.get_positions()[:CHAIN_LENGTH]
final_rg = calculate_radius_of_gyration(chain_positions)
final_contacts = calculate_hydrophobic_contacts(chain_positions, SEQUENCE)

# Calculate statistics
initial_rg = rg_values[0]
mean_rg = np.mean(rg_values[1000:])  # After initial collapse
std_rg = np.std(rg_values[1000:])
mean_contacts = np.mean(contacts_values[1000:])
max_contacts = max(contacts_values)

ke, pe, total = sim.total_energy()

print(f"\nFinal Statistics:")
print(f"  Temperature: {sim.temperature():.1f} K")
print(f"  Initial Rg: {initial_rg:.3f} nm (extended)")
print(f"  Final Rg: {final_rg:.3f} nm (folded)")
print(f"  Collapse ratio: {final_rg/initial_rg:.2f}")
print(f"  Average Rg (equilibrated): {mean_rg:.3f} ± {std_rg:.3f} nm")
print(f"  Hydrophobic contacts: {final_contacts}")
print(f"  Max contacts observed: {max_contacts}")
print(f"  Average contacts: {mean_contacts:.1f}")
print(f"  Total energy: {total:.2f} kJ/mol")

# Calculate compactness
n_hydrophobic = sum(SEQUENCE)
max_possible_contacts = n_hydrophobic * (n_hydrophobic - 1) // 2
print(f"\nFolding Quality:")
print(f"  Contact formation: {final_contacts}/{max_possible_contacts} "
      f"({100*final_contacts/max_possible_contacts:.1f}% of maximum)")

# Save final frame
output_file = os.path.join(os.path.dirname(__file__), '../assets/protein_folding_2d.png')
renderer.save_frame(output_file)
print(f"\nSaved final frame to: {output_file}")

# Create energy summary
from core.renderer import plot_energy_summary

fig = plot_energy_summary(sim, figsize=(12, 8))
summary_file = os.path.join(os.path.dirname(__file__), '../assets/protein_folding_2d_energy.png')
fig.savefig(summary_file, dpi=150, bbox_inches='tight')
print(f"Saved energy summary to: {summary_file}")

# Create folding analysis plot
fig_fold, axes = plt.subplots(2, 2, figsize=(12, 10))

times = np.arange(len(rg_values)) * sim.dt

# Radius of gyration (folding progress)
axes[0, 0].plot(times, rg_values, 'b-', linewidth=1.5)
axes[0, 0].axhline(initial_rg, color='gray', linestyle=':', label='Extended')
axes[0, 0].axhline(mean_rg, color='r', linestyle='--', label=f'Folded: {mean_rg:.3f} nm')
axes[0, 0].set_xlabel('Time (ps)')
axes[0, 0].set_ylabel('Radius of Gyration (nm)')
axes[0, 0].set_title('Folding Progress (Collapse)')
axes[0, 0].legend()
axes[0, 0].grid(True, alpha=0.3)

# Hydrophobic contacts (native structure formation)
axes[0, 1].plot(times, contacts_values, 'r-', linewidth=1.5)
axes[0, 1].axhline(mean_contacts, color='b', linestyle='--',
                   label=f'Average: {mean_contacts:.1f}')
axes[0, 1].set_xlabel('Time (ps)')
axes[0, 1].set_ylabel('Hydrophobic Contacts')
axes[0, 1].set_title('Native Contact Formation')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.3)

# Contact map (final structure)
contact_map = np.zeros((CHAIN_LENGTH, CHAIN_LENGTH))
for i in range(CHAIN_LENGTH):
    for j in range(i+1, CHAIN_LENGTH):
        r = np.linalg.norm(chain_positions[i] - chain_positions[j])
        if r < 0.4:  # Contact cutoff
            contact_map[i, j] = 1
            contact_map[j, i] = 1

im = axes[1, 0].imshow(contact_map, cmap='Blues', interpolation='nearest')
axes[1, 0].set_xlabel('Residue Index')
axes[1, 0].set_ylabel('Residue Index')
axes[1, 0].set_title('Contact Map (Final Structure)')
plt.colorbar(im, ax=axes[1, 0])

# Scatter plot: Contacts vs Rg (folding funnel)
# Color by time to show progression
scatter = axes[1, 1].scatter(rg_values, contacts_values,
                           c=times, cmap='viridis', s=10, alpha=0.6)
axes[1, 1].set_xlabel('Radius of Gyration (nm)')
axes[1, 1].set_ylabel('Hydrophobic Contacts')
axes[1, 1].set_title('Folding Funnel (color = time)')
plt.colorbar(scatter, ax=axes[1, 1], label='Time (ps)')
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
fold_file = os.path.join(os.path.dirname(__file__), '../assets/protein_folding_2d_analysis.png')
fig_fold.savefig(fold_file, dpi=150, bbox_inches='tight')
print(f"Saved folding analysis to: {fold_file}")

print("\nEducational Notes:")
print("  - Hydrophobic residues (red) prefer to be together")
print("  - This drives collapse from extended to compact state")
print("  - Polar residues (blue) are on the surface")
print("  - This is a simplified 2D version of real protein folding")
print("  - Real proteins use 20 amino acid types in 3D!")
print("  - Folding funnel shows path from unfolded to folded state")

renderer.show()
