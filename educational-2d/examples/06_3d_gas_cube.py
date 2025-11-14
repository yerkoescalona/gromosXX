"""
Example 6: 3D Gas in a Cube
============================

Demonstrates:
- 3D molecular dynamics simulation
- Lennard-Jones interactions in 3D
- 3D visualization with matplotlib
- Interactive rotation
- Comparison with 2D case

This shows how the same physics works in 3D,
with 3 degrees of freedom instead of 2.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from core.simulation3d import Simulation3D
from core.renderer3d import Renderer3D, InteractiveRenderer3D

# Simulation parameters
N_PARTICLES = 40
BOX_SIZE = (4.0, 4.0, 4.0)  # nm (cubic box)
TEMPERATURE = 300.0  # K

print("=" * 70)
print("Educational 3D Simulation: Gas in a Cube")
print("=" * 70)
print(f"\nSimulation Parameters:")
print(f"  Number of particles: {N_PARTICLES}")
print(f"  Box size: {BOX_SIZE[0]} x {BOX_SIZE[1]} x {BOX_SIZE[2]} nm")
print(f"  Volume: {BOX_SIZE[0] * BOX_SIZE[1] * BOX_SIZE[2]:.1f} nm³")
print(f"  Temperature: {TEMPERATURE} K")
print(f"  Density: {N_PARTICLES / (BOX_SIZE[0] * BOX_SIZE[1] * BOX_SIZE[2]):.2f} particles/nm³")
print(f"  Periodic boundaries: Yes")
print()

# Create simulation
sim = Simulation3D(
    box_size=BOX_SIZE,
    periodic=True,
    integrator='velocity-verlet',
    dt=0.002
)

# Add particles in a 3D grid
print("Initializing particles in 3D grid...")
grid_size = int(np.ceil(N_PARTICLES**(1/3)))
spacing = min(BOX_SIZE) / grid_size

particle_count = 0
for ix in range(grid_size):
    for iy in range(grid_size):
        for iz in range(grid_size):
            if particle_count >= N_PARTICLES:
                break

            x = (ix + 0.5) * spacing + np.random.uniform(-0.1, 0.1)
            y = (iy + 0.5) * spacing + np.random.uniform(-0.1, 0.1)
            z = (iz + 0.5) * spacing + np.random.uniform(-0.1, 0.1)

            # Wrap to box
            x = x % BOX_SIZE[0]
            y = y % BOX_SIZE[1]
            z = z % BOX_SIZE[2]

            # Different colors based on z-layer for depth perception
            if iz % 3 == 0:
                color = 'blue'
            elif iz % 3 == 1:
                color = 'red'
            else:
                color = 'green'

            sim.add_particle(
                position=[x, y, z],
                mass=1.0,
                atom_type=0,
                color=color
            )
            particle_count += 1

print(f"Added {len(sim.particles)} particles")

# Initialize velocities from Maxwell-Boltzmann distribution (3D)
print(f"Initializing velocities at {TEMPERATURE} K...")
sim.initialize_velocities(TEMPERATURE)

# Enable temperature coupling
sim.set_temperature(TEMPERATURE, tau=0.1)

print(f"Initial temperature: {sim.temperature():.1f} K")

# Create renderer
renderer = Renderer3D(
    box_size=BOX_SIZE,
    figsize=(14, 6),
    azimuth=45,
    elevation=30
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
print("Production Phase with Visualization")
print("=" * 70)
print("\nRotating view during simulation...")
print("Observe particles moving in 3D space!")
print()

for step in range(1000):
    sim.integrate_step()

    # Render every 10 steps with rotating view
    if step % 10 == 0:
        # Slowly rotate the view
        azimuth = 45 + (360 * step / 1000)
        renderer.render(sim, show_bonds=False, azimuth=azimuth)

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

# Theoretical analysis (3D ideal gas)
k_B = 0.00831446261815324  # kJ/(mol·K)
expected_ke_per_particle_3d = 1.5 * k_B * TEMPERATURE  # In 3D: <KE> = 3/2 * k_B * T
actual_ke_per_particle = ke / N_PARTICLES

print(f"\nTheoretical Analysis (3D Ideal Gas):")
print(f"  Expected <KE> per particle: {expected_ke_per_particle_3d:.3f} kJ/mol")
print(f"  Actual <KE> per particle: {actual_ke_per_particle:.3f} kJ/mol")
print(f"  Ratio (should be ~1): {actual_ke_per_particle/expected_ke_per_particle_3d:.3f}")

# Comparison with 2D
expected_ke_per_particle_2d = k_B * TEMPERATURE
print(f"\nComparison with 2D:")
print(f"  2D expected <KE>: {expected_ke_per_particle_2d:.3f} kJ/mol (2 DOF)")
print(f"  3D expected <KE>: {expected_ke_per_particle_3d:.3f} kJ/mol (3 DOF)")
print(f"  Ratio 3D/2D: {expected_ke_per_particle_3d/expected_ke_per_particle_2d:.2f}")

# Save final frame
output_file = os.path.join(os.path.dirname(__file__), '../assets/gas_3d_cube.png')
renderer.save_frame(output_file)
print(f"\nSaved final frame to: {output_file}")

# Create a rotating movie
print("\n" + "=" * 70)
print("Creating Rotating 3D Movie")
print("=" * 70)

try:
    from core.movie_utils import make_rotating_3d_movie

    # Take a snapshot and create rotating view
    output_movie = os.path.join(os.path.dirname(__file__), '../assets/gas_3d_rotation.mp4')

    print("\nCreating 360° rotation movie...")
    make_rotating_3d_movie(
        sim, renderer,
        n_frames=180,  # Half rotation
        output_file=output_movie,
        elevation=30,
        fps=30
    )
    print(f"Saved rotating movie to: {output_movie}")

except Exception as e:
    print(f"Could not create rotation movie: {e}")
    print("Install ffmpeg for video export")

# Optional: Interactive 3D view
print("\n" + "=" * 70)
print("Interactive 3D View")
print("=" * 70)
print("\nLaunching interactive 3D viewer...")
print("You can:")
print("  - Click and drag to rotate")
print("  - Scroll to zoom")
print("  - Close window to continue")
print()

try:
    interactive_renderer = InteractiveRenderer3D(box_size=BOX_SIZE, figsize=(10, 8))

    # Run a few more steps with interactive view
    for step in range(200):
        sim.integrate_step()

        if step % 5 == 0:
            interactive_renderer.render(sim, show_bonds=False)

    interactive_renderer.show()

except Exception as e:
    print(f"Could not create interactive view: {e}")

print("\n" + "=" * 70)
print("Educational Notes:")
print("=" * 70)
print("""
3D vs 2D Molecular Dynamics:

1. Degrees of Freedom:
   - 2D: Each particle has 2 DOF (x, y)
   - 3D: Each particle has 3 DOF (x, y, z)

2. Equipartition Theorem:
   - 2D: <KE> = k_B * T (total for 2 DOF)
   - 3D: <KE> = 3/2 * k_B * T (total for 3 DOF)
   - Per DOF: <KE_per_DOF> = 1/2 * k_B * T (same in 2D and 3D!)

3. Spatial Distribution:
   - 2D: Particles confined to a plane
   - 3D: Particles fill a volume

4. Collisions:
   - 3D has more "room" - particles collide less frequently
   - Same density in 3D = particles are more spread out

5. Visualization:
   - 2D: Easy to see all particles at once
   - 3D: Need rotation to see full structure
   - Depth perception via colors helps in 3D

This demonstrates that the same physics (LJ potential, integration,
thermostat) works identically in 2D and 3D - only the geometry changes!
""")

renderer.show()
