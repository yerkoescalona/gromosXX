"""
Example 3: Crystal Lattice Formation
=====================================

Demonstrates:
- Hexagonal close-packed lattice (2D equivalent of FCC)
- Lattice vibrations (phonons)
- Melting transition
- Order-disorder transition
- Crystalline vs liquid vs gas phases

This simulation shows how particles arrange themselves in ordered
structures at low temperatures, and how thermal energy disrupts
this order, leading to melting.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from core import Simulation2D, Renderer2D

# Simulation parameters
LATTICE_SIZE = 8  # 8x8 lattice
N_PARTICLES = LATTICE_SIZE * LATTICE_SIZE
BOX_SIZE = (5.0, 5.0)
INITIAL_TEMPERATURE = 50.0  # K (low for crystal)

print("=" * 70)
print("Educational 2D Simulation: Crystal Lattice")
print("=" * 70)
print(f"\nSimulation Parameters:")
print(f"  Lattice size: {LATTICE_SIZE} x {LATTICE_SIZE}")
print(f"  Number of particles: {N_PARTICLES}")
print(f"  Box size: {BOX_SIZE[0]} x {BOX_SIZE[1]} nm")
print(f"  Initial temperature: {INITIAL_TEMPERATURE} K")
print(f"  Periodic boundaries: Yes")
print()

# Create simulation
sim = Simulation2D(
    box_size=BOX_SIZE,
    periodic=True,
    integrator='velocity-verlet',
    dt=0.001  # Smaller timestep for stability
)

# Adjust LJ parameters for crystal formation
# Balance between attraction and repulsion
sim.lj_params[0] = {
    'C6': 0.003,
    'C12': 0.00003,
    'sigma': 0.35,
    'epsilon': 0.8
}

print("Initializing particles in perfect hexagonal lattice...")

# Create hexagonal close-packed lattice
# This is the 2D equivalent of FCC (most efficient packing)
spacing = BOX_SIZE[0] / LATTICE_SIZE
particle_idx = 0

for i in range(LATTICE_SIZE):
    for j in range(LATTICE_SIZE):
        # Hexagonal lattice with offset every other row
        x = i * spacing + 0.5 * spacing
        y = j * spacing + 0.5 * spacing

        # Hexagonal offset
        if j % 2 == 1:
            x += 0.5 * spacing

        # Wrap to periodic box
        x = x % BOX_SIZE[0]
        y = y % BOX_SIZE[1]

        # Color based on position for visual tracking
        if (i + j) % 2 == 0:
            color = 'blue'
        else:
            color = 'lightblue'

        sim.add_particle(
            position=[x, y],
            mass=1.0,
            atom_type=0,
            color=color
        )
        particle_idx += 1

print(f"Created {len(sim.particles)} particles in lattice")

# Initialize velocities at low temperature
print(f"Initializing velocities at {INITIAL_TEMPERATURE} K...")
sim.initialize_velocities(INITIAL_TEMPERATURE)

# Enable temperature coupling
sim.set_temperature(INITIAL_TEMPERATURE, tau=0.1)

print(f"Initial temperature: {sim.temperature():.1f} K")

# Create renderer
renderer = Renderer2D(
    box_size=BOX_SIZE,
    figsize=(14, 6),
    show_forces=False,
    show_velocities=False,
    show_trails=False
)

# Phase 1: Equilibrate crystal at low temperature
print("\n" + "=" * 70)
print("Phase 1: Crystal Equilibration (1000 steps at 50 K)")
print("=" * 70)
print("\nObserve: Small vibrations around lattice sites (phonons)")
print()

for step in range(1000):
    sim.integrate_step()

    if step % 10 == 0:
        renderer.render(sim, show_bonds=False)

    if step % 200 == 0:
        ke, pe, total = sim.total_energy()
        print(f"Step {step:4d}: T = {sim.temperature():6.1f} K, "
              f"E_tot = {total:8.2f} kJ/mol")

# Save crystal state
output_crystal = os.path.join(os.path.dirname(__file__), '../assets/crystal_lattice_ordered.png')
renderer.save_frame(output_crystal)
print(f"\nSaved crystal state to: {output_crystal}")

# Phase 2: Gradual heating (melting transition)
print("\n" + "=" * 70)
print("Phase 2: Heating and Melting (2000 steps)")
print("=" * 70)
print("\nRamping temperature from 50 K to 300 K...")
print("Observe: Lattice breaks down as particles gain thermal energy")
print()

initial_positions = sim.get_positions().copy()

for step in range(2000):
    # Gradually increase temperature
    current_temp = 50.0 + (250.0 * step / 2000)
    sim.set_temperature(current_temp, tau=0.1)

    sim.integrate_step()

    if step % 10 == 0:
        renderer.render(sim, show_bonds=False)

    if step % 200 == 0:
        ke, pe, total = sim.total_energy()
        current_positions = sim.get_positions()

        # Calculate mean square displacement from initial lattice
        msd = np.mean(np.sum((current_positions - initial_positions)**2, axis=1))

        print(f"Step {step:4d}: T = {sim.temperature():6.1f} K, "
              f"E_tot = {total:8.2f} kJ/mol, MSD = {msd:.4f} nm²")

# Save melted state
output_melted = os.path.join(os.path.dirname(__file__), '../assets/crystal_lattice_melted.png')
renderer.save_frame(output_melted)
print(f"\nSaved melted state to: {output_melted}")

# Phase 3: High temperature (liquid/gas phase)
print("\n" + "=" * 70)
print("Phase 3: High Temperature Dynamics (1000 steps at 300 K)")
print("=" * 70)
print("\nObserve: Particles freely diffuse (liquid-like or gas-like)")
print()

for step in range(1000):
    sim.integrate_step()

    if step % 10 == 0:
        renderer.render(sim, show_bonds=False)

    if step % 200 == 0:
        ke, pe, total = sim.total_energy()
        print(f"Step {step:4d}: T = {sim.temperature():6.1f} K, "
              f"E_tot = {total:8.2f} kJ/mol")

print("\n" + "=" * 70)
print("Simulation Complete!")
print("=" * 70)

# Final analysis
ke, pe, total = sim.total_energy()
final_temp = sim.temperature()
final_positions = sim.get_positions()

print(f"\nFinal Statistics:")
print(f"  Temperature: {final_temp:.1f} K")
print(f"  Kinetic Energy: {ke:.2f} kJ/mol")
print(f"  Potential Energy: {pe:.2f} kJ/mol")
print(f"  Total Energy: {total:.2f} kJ/mol")

# Calculate order parameter (deviation from initial lattice)
msd_final = np.mean(np.sum((final_positions - initial_positions)**2, axis=1))
print(f"  Mean square displacement from lattice: {msd_final:.4f} nm²")

# Save final frame
output_file = os.path.join(os.path.dirname(__file__), '../assets/crystal_lattice_final.png')
renderer.save_frame(output_file)
print(f"\nSaved final frame to: {output_file}")

# Create energy summary
from core.renderer import plot_energy_summary

fig = plot_energy_summary(sim, figsize=(12, 8))
summary_file = os.path.join(os.path.dirname(__file__), '../assets/crystal_lattice_energy.png')
fig.savefig(summary_file, dpi=150, bbox_inches='tight')
print(f"Saved energy summary to: {summary_file}")

# Create order parameter plot
fig_order, ax = plt.subplots(figsize=(10, 6))

# Calculate MSD over time (if we had stored positions)
# For now, show energy components as proxy for order
times = np.arange(len(sim.energies['total'])) * sim.dt

ax.plot(times, sim.energies['potential'], 'r-', linewidth=2, label='Potential Energy')
ax.axvline(1.0, color='k', linestyle='--', alpha=0.5, label='Start Heating')
ax.axvline(3.0, color='k', linestyle='--', alpha=0.5, label='End Heating')
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Potential Energy (kJ/mol)')
ax.set_title('Melting Transition (PE decreases as lattice breaks)')
ax.legend()
ax.grid(True, alpha=0.3)

order_file = os.path.join(os.path.dirname(__file__), '../assets/crystal_lattice_order.png')
fig_order.savefig(order_file, dpi=150, bbox_inches='tight')
print(f"Saved order parameter plot to: {order_file}")

print("\nEducational Notes:")
print("  - At low T: Particles vibrate around lattice sites (solid)")
print("  - During heating: Lattice breaks down (melting)")
print("  - At high T: Particles diffuse freely (liquid/gas)")
print("  - This demonstrates first-order phase transition!")

renderer.show()
