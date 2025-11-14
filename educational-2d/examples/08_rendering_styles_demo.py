"""
Example 8: Rendering Styles Demo
=================================

Demonstrates different ways to visualize particles/molecules:
- Multiple rendering styles (spheres, halos, wireframes, etc.)
- Color schemes (velocity, energy, temperature, etc.)
- Size scaling options
- Ball-and-stick vs space-filling models

Shows all the visual representations available for educational simulations.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from core import Simulation2D, Renderer2D
from core.render_styles import (
    ParticleRenderer, RenderStyle, ColorScheme, create_colorbar
)

print("=" * 70)
print("Educational Simulation: Rendering Styles Demo")
print("=" * 70)

# Create a simple simulation
N_PARTICLES = 30
BOX_SIZE = (4.0, 4.0)
TEMPERATURE = 300.0

print(f"\nCreating simulation with {N_PARTICLES} particles...")

sim = Simulation2D(
    box_size=BOX_SIZE,
    periodic=True,
    integrator='velocity-verlet',
    dt=0.002
)

# Add particles with variety of velocities
grid_size = int(np.ceil(np.sqrt(N_PARTICLES)))
spacing = min(BOX_SIZE) / grid_size

for i in range(N_PARTICLES):
    ix = i % grid_size
    iy = i // grid_size
    x = (ix + 0.5) * spacing
    y = (iy + 0.5) * spacing

    # Assign different atom types
    atom_type = i % 3
    colors = ['blue', 'red', 'green']

    sim.add_particle(
        position=[x, y],
        mass=1.0,
        atom_type=atom_type,
        color=colors[atom_type]
    )

# Add some bonds to show ball-and-stick style
print("Adding bonds for molecular structure...")
for i in range(N_PARTICLES - 1):
    if i % 5 != 4:  # Create chains of 5 particles
        sim.add_bond(i=i, j=i+1, k=1000.0, r0=0.5)

# Initialize with varied velocities
sim.initialize_velocities(TEMPERATURE)
sim.set_temperature(TEMPERATURE, tau=0.1)

# Run brief equilibration
print("Equilibrating...")
for step in range(200):
    sim.integrate_step()

print("\n" + "=" * 70)
print("Demonstrating Different Rendering Styles")
print("=" * 70)

# Create figure with multiple subplots
fig = plt.figure(figsize=(18, 12))

styles_to_demo = [
    (RenderStyle.SPHERES, ColorScheme.ATOM_TYPE, "Standard Spheres (by atom type)"),
    (RenderStyle.HALOS, ColorScheme.ATOM_TYPE, "Halos (glowing spheres)"),
    (RenderStyle.WIREFRAME, ColorScheme.ATOM_TYPE, "Wireframe Circles"),
    (RenderStyle.VELOCITY_ARROWS, ColorScheme.ATOM_TYPE, "Velocity Arrows"),
    (RenderStyle.SPHERES, ColorScheme.VELOCITY, "Color by Velocity"),
    (RenderStyle.SPHERES, ColorScheme.KINETIC_ENERGY, "Color by Kinetic Energy"),
    (RenderStyle.SPHERES, ColorScheme.TEMPERATURE, "Color by Temperature"),
    (RenderStyle.BALLS_AND_STICKS, ColorScheme.ATOM_TYPE, "Ball-and-Stick Model"),
    (RenderStyle.SPACE_FILLING, ColorScheme.ATOM_TYPE, "Space-Filling Model"),
]

for idx, (style, color_scheme, title) in enumerate(styles_to_demo):
    ax = fig.add_subplot(3, 3, idx + 1)

    print(f"\nRendering: {title}")

    # Create renderer with specific style
    renderer = ParticleRenderer(
        style=style,
        color_scheme=color_scheme,
        size_scale=100.0,
        alpha=0.8,
        colormap='plasma'
    )

    # Setup axes
    ax.set_xlim(0, BOX_SIZE[0])
    ax.set_ylim(0, BOX_SIZE[1])
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.set_xlabel('x (nm)', fontsize=8)
    ax.set_ylabel('y (nm)', fontsize=8)

    # Draw box
    from matplotlib.patches import Rectangle
    rect = Rectangle((0, 0), BOX_SIZE[0], BOX_SIZE[1],
                     linewidth=2, edgecolor='black', facecolor='none')
    ax.add_patch(rect)

    # Render bonds first (if applicable)
    if style in [RenderStyle.BALLS_AND_STICKS, RenderStyle.SPACE_FILLING]:
        renderer.render_bonds(ax, sim, bond_thickness=3.0 if style == RenderStyle.SPACE_FILLING else 2.0)

    # Render particles
    renderer.render_particles_2d(ax, sim)

    # Add colorbar for continuous color schemes
    if color_scheme in [ColorScheme.VELOCITY, ColorScheme.KINETIC_ENERGY,
                        ColorScheme.TEMPERATURE]:
        create_colorbar(ax, renderer, sim)

    # Add info text
    info_text = f"Particles: {len(sim.particles)}\n"
    info_text += f"Bonds: {len(sim.bonds)}\n"
    info_text += f"T: {sim.temperature():.1f} K"

    ax.text(0.02, 0.98, info_text,
           transform=ax.transAxes,
           verticalalignment='top',
           fontsize=7,
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()

# Save the comparison figure
output_file = os.path.join(os.path.dirname(__file__), '../assets/rendering_styles_comparison.png')
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"\nSaved rendering styles comparison to: {output_file}")

plt.show()

# Now demonstrate animated styles
print("\n" + "=" * 70)
print("Animated Demonstration: Velocity-Colored Particles")
print("=" * 70)

fig2, ax2 = plt.subplots(figsize=(10, 10))

# Use velocity coloring
renderer_anim = ParticleRenderer(
    style=RenderStyle.HALOS,
    color_scheme=ColorScheme.VELOCITY,
    size_scale=150.0,
    alpha=0.8,
    colormap='hot'
)

print("\nRunning animated simulation with velocity coloring...")
print("Colors show particle speeds (blue=slow, red=fast)")

for step in range(500):
    sim.integrate_step()

    if step % 10 == 0:
        ax2.clear()

        # Setup axes
        ax2.set_xlim(0, BOX_SIZE[0])
        ax2.set_ylim(0, BOX_SIZE[1])
        ax2.set_aspect('equal')
        ax2.set_title('Velocity-Colored Halos (speeds shown by color)', fontsize=14)
        ax2.set_xlabel('x (nm)')
        ax2.set_ylabel('y (nm)')

        # Draw box
        rect = Rectangle((0, 0), BOX_SIZE[0], BOX_SIZE[1],
                        linewidth=2, edgecolor='black', facecolor='none')
        ax2.add_patch(rect)

        # Render
        renderer_anim.render_particles_2d(ax2, sim)
        create_colorbar(ax2, renderer_anim, sim, label="Speed")

        # Info
        info = f"Step: {step}\nT: {sim.temperature():.1f} K\n"
        velocities = sim.get_velocities()
        speeds = np.linalg.norm(velocities, axis=1)
        info += f"Max speed: {speeds.max():.3f} nm/ps\n"
        info += f"Avg speed: {speeds.mean():.3f} nm/ps"

        ax2.text(0.02, 0.98, info,
                transform=ax2.transAxes,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7),
                fontsize=10)

        plt.draw()
        plt.pause(0.001)

print("\n" + "=" * 70)
print("Rendering Styles Summary")
print("=" * 70)

print("""
Available Rendering Styles:
---------------------------
1. SPHERES - Standard filled circles (default)
2. HALOS - Spheres with glowing outer ring
3. WIREFRAME - Hollow circles (outline only)
4. VELOCITY_ARROWS - Particles with velocity vectors
5. BALLS_AND_STICKS - Small spheres + bonds (molecular)
6. SPACE_FILLING - Large touching spheres (CPK model)
7. VDW_RADII - Sized by Van der Waals radii
8. POINTS - Tiny dots (for many particles)

Available Color Schemes:
-----------------------
1. ATOM_TYPE - Fixed colors by atom type (default)
2. VELOCITY - Color by speed (blue=slow, red=fast)
3. KINETIC_ENERGY - Color by kinetic energy
4. TEMPERATURE - Color by local temperature
5. FORCE - Color by force magnitude
6. CPK - Standard CPK molecular colors
7. RAINBOW - Rainbow gradient by index
8. CUSTOM - User-provided colors

Usage Examples:
--------------
# Halos with velocity coloring
renderer = ParticleRenderer(
    style=RenderStyle.HALOS,
    color_scheme=ColorScheme.VELOCITY,
    colormap='plasma'
)

# Ball-and-stick model
renderer = ParticleRenderer(
    style=RenderStyle.BALLS_AND_STICKS,
    color_scheme=ColorScheme.ATOM_TYPE,
    size_scale=80
)

# Temperature visualization
renderer = ParticleRenderer(
    style=RenderStyle.SPHERES,
    color_scheme=ColorScheme.TEMPERATURE,
    colormap='coolwarm',
    alpha=0.9
)

Integration with Existing Renderers:
-----------------------------------
These styles can be used with Renderer2D, Renderer3D, or standalone.
Perfect for creating publication-quality figures with different
visual representations of the same molecular system!
""")

plt.show()
