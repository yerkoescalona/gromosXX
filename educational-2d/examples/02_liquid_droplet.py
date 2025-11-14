"""
Example 2: Liquid Droplet Formation
====================================

Demonstrates:
- Cohesive forces (attractive LJ interactions)
- Surface tension
- Droplet coalescence
- Phase separation
- No periodic boundaries (vacuum)

This simulation shows how particles with attractive interactions
spontaneously form liquid-like droplets in vacuum, demonstrating
the balance between cohesive forces and kinetic energy.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from core import Simulation2D, Renderer2D

# Simulation parameters
N_PARTICLES = 80
BOX_SIZE = (8.0, 8.0)  # nm (larger box for droplet formation)
TEMPERATURE = 100.0    # K (low temperature for condensation)

print("=" * 70)
print("Educational 2D Simulation: Liquid Droplet Formation")
print("=" * 70)
print(f"\nSimulation Parameters:")
print(f"  Number of particles: {N_PARTICLES}")
print(f"  Box size: {BOX_SIZE[0]} x {BOX_SIZE[1]} nm")
print(f"  Temperature: {TEMPERATURE} K (low for condensation)")
print(f"  Periodic boundaries: No (vacuum)")
print(f"  Attractive interactions: Strong (liquid-like)")
print()

# Create simulation (no periodic boundaries for droplet in vacuum)
sim = Simulation2D(
    box_size=BOX_SIZE,
    periodic=False,
    integrator='velocity-verlet',
    dt=0.002
)

# Modify LJ parameters for stronger attraction (liquid-like)
# Increase epsilon (well depth) for more cohesion
sim.lj_params[0] = {
    'C6': 0.005,    # Increased attraction
    'C12': 0.00001,
    'sigma': 0.3,
    'epsilon': 1.5   # Stronger interactions
}

print("Initializing particles in two separate clusters...")

# Create two small clusters that will coalesce
cluster1_center = [2.5, 4.0]
cluster2_center = [5.5, 4.0]

cluster_size = N_PARTICLES // 2

for i in range(N_PARTICLES):
    if i < cluster_size:
        # First cluster (blue)
        angle = 2 * np.pi * i / cluster_size
        radius = 0.8 * np.random.uniform(0.3, 0.8)
        x = cluster1_center[0] + radius * np.cos(angle)
        y = cluster1_center[1] + radius * np.sin(angle)
        color = 'blue'
    else:
        # Second cluster (red)
        angle = 2 * np.pi * (i - cluster_size) / cluster_size
        radius = 0.8 * np.random.uniform(0.3, 0.8)
        x = cluster2_center[0] + radius * np.cos(angle)
        y = cluster2_center[1] + radius * np.sin(angle)
        color = 'red'

    sim.add_particle(
        position=[x, y],
        mass=1.0,
        atom_type=0,
        color=color
    )

# Initialize velocities at low temperature
print(f"Initializing velocities at {TEMPERATURE} K...")
sim.initialize_velocities(TEMPERATURE)

# Light temperature coupling to remove excess energy
sim.set_temperature(TEMPERATURE, tau=0.5)

print(f"Initial temperature: {sim.temperature():.1f} K")

# Create renderer with trails to show motion
renderer = Renderer2D(
    box_size=BOX_SIZE,
    figsize=(14, 6),
    show_forces=False,
    show_velocities=False,
    show_trails=True,
    trail_length=30
)

# Run simulation with visualization
print("\n" + "=" * 70)
print("Watching Droplet Coalescence (3000 steps)")
print("=" * 70)
print("\nObserve how:")
print("  1. Particles attract each other (LJ potential)")
print("  2. Two clusters move toward each other")
print("  3. Clusters merge into a single droplet")
print("  4. Surface tension creates a round shape")
print()

for step in range(3000):
    sim.integrate_step()

    # Render every 5 steps
    if step % 5 == 0:
        renderer.render(sim, show_bonds=False)

    if step % 300 == 0:
        ke, pe, total = sim.total_energy()
        print(f"Step {step:4d}: T = {sim.temperature():6.1f} K, "
              f"E_pot = {pe:8.2f} kJ/mol, E_tot = {total:8.2f} kJ/mol")

        # Calculate radius of gyration (measure of cluster size)
        positions = sim.get_positions()
        com = np.mean(positions, axis=0)
        distances = np.linalg.norm(positions - com, axis=1)
        rg = np.sqrt(np.mean(distances**2))
        print(f"         Radius of gyration: {rg:.3f} nm")

print("\n" + "=" * 70)
print("Simulation Complete!")
print("=" * 70)

# Final analysis
positions = sim.get_positions()
com = np.mean(positions, axis=0)
distances = np.linalg.norm(positions - com, axis=1)
rg = np.sqrt(np.mean(distances**2))
max_dist = np.max(distances)

ke, pe, total = sim.total_energy()
final_temp = sim.temperature()

print(f"\nFinal Statistics:")
print(f"  Temperature: {final_temp:.1f} K")
print(f"  Potential Energy: {pe:.2f} kJ/mol")
print(f"  Center of mass: ({com[0]:.2f}, {com[1]:.2f}) nm")
print(f"  Radius of gyration: {rg:.3f} nm")
print(f"  Maximum distance from COM: {max_dist:.3f} nm")
print(f"  Cohesive energy per particle: {pe/N_PARTICLES:.3f} kJ/mol")

# Save final frame
output_file = os.path.join(os.path.dirname(__file__), '../assets/liquid_droplet.png')
renderer.save_frame(output_file)
print(f"\nSaved final frame to: {output_file}")

# Create energy summary
from core.renderer import plot_energy_summary

fig = plot_energy_summary(sim, figsize=(12, 8))
summary_file = os.path.join(os.path.dirname(__file__), '../assets/liquid_droplet_energy.png')
fig.savefig(summary_file, dpi=150, bbox_inches='tight')
print(f"Saved energy summary to: {summary_file}")

# Additional analysis plot: radial distribution
fig_rdf, ax = plt.subplots(figsize=(10, 6))

# Calculate pairwise distances
n = len(positions)
distances_all = []
for i in range(n):
    for j in range(i+1, n):
        r = np.linalg.norm(positions[i] - positions[j])
        distances_all.append(r)

# Histogram of distances
hist, bins = np.histogram(distances_all, bins=50, range=(0, 3))
bin_centers = (bins[:-1] + bins[1:]) / 2

ax.bar(bin_centers, hist, width=bins[1]-bins[0], alpha=0.7, edgecolor='black')
ax.set_xlabel('Distance (nm)')
ax.set_ylabel('Count')
ax.set_title('Pair Distance Distribution (Radial Distribution Function)')
ax.grid(True, alpha=0.3)

rdf_file = os.path.join(os.path.dirname(__file__), '../assets/liquid_droplet_rdf.png')
fig_rdf.savefig(rdf_file, dpi=150, bbox_inches='tight')
print(f"Saved RDF plot to: {rdf_file}")

renderer.show()
