"""
Example 1: Gas in a Box
========================

Demonstrates:
- Lennard-Jones interactions
- Periodic boundary conditions
- Maxwell-Boltzmann velocity distribution
- Ideal gas behavior in 2D
- Energy conservation

This simulation shows how particles with purely repulsive/attractive
interactions (no bonds) behave in a confined space, similar to
a simple gas or liquid depending on temperature and density.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Add educational-2d to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from core import Simulation2D, Renderer2D

# Simulation parameters
N_PARTICLES = 50
BOX_SIZE = (5.0, 5.0)  # nm
TEMPERATURE = 300.0    # K
DENSITY = N_PARTICLES / (BOX_SIZE[0] * BOX_SIZE[1])  # particles/nm²

print("=" * 70)
print("Educational 2D Simulation: Gas in a Box")
print("=" * 70)
print(f"\nSimulation Parameters:")
print(f"  Number of particles: {N_PARTICLES}")
print(f"  Box size: {BOX_SIZE[0]} x {BOX_SIZE[1]} nm")
print(f"  Temperature: {TEMPERATURE} K")
print(f"  Density: {DENSITY:.2f} particles/nm²")
print(f"  Periodic boundaries: Yes")
print()

# Create simulation
sim = Simulation2D(
    box_size=BOX_SIZE,
    periodic=True,
    integrator='velocity-verlet',
    dt=0.002  # 2 fs timestep
)

# Add particles in a grid with random perturbations
grid_size = int(np.ceil(np.sqrt(N_PARTICLES)))
spacing = min(BOX_SIZE) / grid_size

print("Initializing particles...")
for i in range(N_PARTICLES):
    ix = i % grid_size
    iy = i // grid_size

    # Grid position with random perturbation
    x = (ix + 0.5) * spacing + np.random.uniform(-0.1, 0.1)
    y = (iy + 0.5) * spacing + np.random.uniform(-0.1, 0.1)

    # Wrap to box
    x = x % BOX_SIZE[0]
    y = y % BOX_SIZE[1]

    # Different colors for visual interest
    if i % 3 == 0:
        color = 'blue'
        atom_type = 0
    elif i % 3 == 1:
        color = 'red'
        atom_type = 1
    else:
        color = 'green'
        atom_type = 2

    sim.add_particle(
        position=[x, y],
        mass=1.0,
        atom_type=atom_type,
        color=color
    )

# Initialize velocities from Maxwell-Boltzmann distribution
print(f"Initializing velocities at {TEMPERATURE} K...")
sim.initialize_velocities(TEMPERATURE)

# Enable temperature coupling (Berendsen thermostat)
sim.set_temperature(TEMPERATURE, tau=0.1)

print(f"Initial temperature: {sim.temperature():.1f} K")

# Create renderer
renderer = Renderer2D(
    box_size=BOX_SIZE,
    figsize=(14, 6),
    show_forces=False,
    show_velocities=False,
    show_trails=False
)

# Equilibration phase
print("\n" + "=" * 70)
print("Equilibration Phase (1000 steps)")
print("=" * 70)

for step in range(1000):
    sim.integrate_step()

    if step % 200 == 0:
        print(f"Step {step:4d}: T = {sim.temperature():6.1f} K, "
              f"E_total = {sim.total_energy()[2]:8.2f} kJ/mol")

# Production phase with visualization
print("\n" + "=" * 70)
print("Production Phase (2000 steps with visualization)")
print("=" * 70)

for step in range(2000):
    sim.integrate_step()

    # Render every 10 steps
    if step % 10 == 0:
        renderer.render(sim, show_bonds=False)

    if step % 200 == 0:
        print(f"Step {step:4d}: T = {sim.temperature():6.1f} K, "
              f"E_total = {sim.total_energy()[2]:8.2f} kJ/mol")

print("\n" + "=" * 70)
print("Simulation Complete!")
print("=" * 70)

# Final statistics
ke, pe, total = sim.total_energy()
final_temp = sim.temperature()

print(f"\nFinal Statistics:")
print(f"  Temperature: {final_temp:.1f} K")
print(f"  Kinetic Energy: {ke:.2f} kJ/mol")
print(f"  Potential Energy: {pe:.2f} kJ/mol")
print(f"  Total Energy: {total:.2f} kJ/mol")
print(f"  Average LJ energy per particle: {pe/N_PARTICLES:.3f} kJ/mol")

# Calculate average kinetic energy per particle (should be ~k_B*T in 2D)
k_B = 0.00831446261815324  # kJ/(mol·K)
expected_ke_per_particle = k_B * TEMPERATURE  # In 2D: <KE> = k_B*T per particle
actual_ke_per_particle = ke / N_PARTICLES

print(f"\nTheoretical Analysis (2D Ideal Gas):")
print(f"  Expected <KE> per particle: {expected_ke_per_particle:.3f} kJ/mol")
print(f"  Actual <KE> per particle: {actual_ke_per_particle:.3f} kJ/mol")
print(f"  Ratio (should be ~1): {actual_ke_per_particle/expected_ke_per_particle:.3f}")

# Save final frame
output_file = os.path.join(os.path.dirname(__file__), '../assets/gas_in_box.png')
renderer.save_frame(output_file)
print(f"\nSaved final frame to: {output_file}")

# Create energy summary plot
from core.renderer import plot_energy_summary

fig = plot_energy_summary(sim, figsize=(12, 8))
summary_file = os.path.join(os.path.dirname(__file__), '../assets/gas_in_box_energy.png')
fig.savefig(summary_file, dpi=150, bbox_inches='tight')
print(f"Saved energy summary to: {summary_file}")

# Show plots
renderer.show()
