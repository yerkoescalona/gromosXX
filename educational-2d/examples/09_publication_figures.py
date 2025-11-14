"""
Example 9: Creating Publication-Quality Figures
================================================

Demonstrates how to create high-quality figures for publications,
presentations, and educational materials using different rendering styles.

Features:
- High-resolution output
- Multiple panel figures
- Professional color schemes
- Annotations and labels
- Comparison views
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from core import Simulation2D
from core.render_styles import (
    ParticleRenderer, RenderStyle, ColorScheme, create_colorbar
)
from matplotlib.patches import Rectangle, FancyBboxPatch

print("=" * 70)
print("Creating Publication-Quality Figures")
print("=" * 70)

# Set publication-quality matplotlib parameters
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'sans-serif',
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 150,
})

# Create three different systems for comparison
print("\nCreating three different molecular systems...")

# System 1: Low temperature (ordered)
sim_cold = Simulation2D(box_size=(3.0, 3.0), periodic=True, dt=0.002)
for i in range(25):
    ix, iy = i % 5, i // 5
    x, y = (ix + 0.5) * 0.6, (iy + 0.5) * 0.6
    sim_cold.add_particle(position=[x, y], color='blue', atom_type=0)
sim_cold.initialize_velocities(50.0)
sim_cold.set_temperature(50.0, tau=0.1)

# System 2: Medium temperature
sim_warm = Simulation2D(box_size=(3.0, 3.0), periodic=True, dt=0.002)
for i in range(25):
    ix, iy = i % 5, i // 5
    x, y = (ix + 0.5) * 0.6, (iy + 0.5) * 0.6
    sim_warm.add_particle(position=[x, y], color='orange', atom_type=1)
sim_warm.initialize_velocities(200.0)
sim_warm.set_temperature(200.0, tau=0.1)

# System 3: High temperature (disordered)
sim_hot = Simulation2D(box_size=(3.0, 3.0), periodic=True, dt=0.002)
for i in range(25):
    ix, iy = i % 5, i // 5
    x, y = (ix + 0.5) * 0.6, (iy + 0.5) * 0.6
    sim_hot.add_particle(position=[x, y], color='red', atom_type=2)
sim_hot.initialize_velocities(500.0)
sim_hot.set_temperature(500.0, tau=0.1)

# Equilibrate
print("Equilibrating systems...")
for sim in [sim_cold, sim_warm, sim_hot]:
    for _ in range(1000):
        sim.integrate_step()

print("Systems equilibrated!")

# ===== FIGURE 1: Temperature Comparison =====
print("\n" + "=" * 70)
print("Figure 1: Temperature Effects on Molecular Motion")
print("=" * 70)

fig1 = plt.figure(figsize=(16, 5))
gs = GridSpec(1, 3, figure=fig1, wspace=0.3)

systems = [
    (sim_cold, "Low Temperature (50 K)", RenderStyle.SPHERES),
    (sim_warm, "Medium Temperature (200 K)", RenderStyle.SPHERES),
    (sim_hot, "High Temperature (500 K)", RenderStyle.SPHERES)
]

for idx, (sim, title, style) in enumerate(systems):
    ax = fig1.add_subplot(gs[0, idx])

    # Use velocity coloring
    renderer = ParticleRenderer(
        style=style,
        color_scheme=ColorScheme.VELOCITY,
        size_scale=150.0,
        alpha=0.85,
        colormap='plasma'
    )

    # Setup
    ax.set_xlim(0, 3.0)
    ax.set_ylim(0, 3.0)
    ax.set_aspect('equal')
    ax.set_title(title, fontweight='bold', pad=10)
    ax.set_xlabel('Position x (nm)')
    if idx == 0:
        ax.set_ylabel('Position y (nm)')

    # Draw box
    rect = Rectangle((0, 0), 3.0, 3.0, linewidth=2.5,
                     edgecolor='black', facecolor='none')
    ax.add_patch(rect)

    # Render
    renderer.render_particles_2d(ax, sim)

    # Add colorbar
    create_colorbar(ax, renderer, sim, label="Speed (nm/ps)")

    # Statistics box
    velocities = sim.get_velocities()
    speeds = np.linalg.norm(velocities, axis=1)

    stats_text = f"T = {sim.temperature():.0f} K\n"
    stats_text += f"⟨v⟩ = {speeds.mean():.3f} nm/ps\n"
    stats_text += f"v_max = {speeds.max():.3f} nm/ps"

    ax.text(0.05, 0.95, stats_text,
           transform=ax.transAxes,
           verticalalignment='top',
           fontsize=9,
           bbox=dict(boxstyle='round', facecolor='white',
                    edgecolor='black', alpha=0.8, linewidth=1.5))

    # Add panel label
    ax.text(0.02, 0.02, f'({chr(97+idx)})',
           transform=ax.transAxes,
           fontsize=14, fontweight='bold',
           bbox=dict(boxstyle='circle', facecolor='yellow',
                    edgecolor='black', alpha=0.8))

fig1.suptitle('Effect of Temperature on Molecular Velocity Distribution',
             fontsize=15, fontweight='bold', y=0.98)

output1 = os.path.join(os.path.dirname(__file__), '../assets/pub_fig1_temperature.png')
fig1.savefig(output1, dpi=300, bbox_inches='tight')
print(f"Saved Figure 1 to: {output1}")

# ===== FIGURE 2: Rendering Style Comparison =====
print("\n" + "=" * 70)
print("Figure 2: Visual Representation Styles")
print("=" * 70)

# Create a small molecular system
sim_mol = Simulation2D(box_size=(3.0, 3.0), periodic=False, dt=0.0005)

# Create a water-like molecule
center = np.array([1.5, 1.5])
sim_mol.add_particle(position=center, color='red', mass=16.0)  # O
sim_mol.add_particle(position=center + [0.15, 0.1], color='lightblue', mass=1.0)  # H1
sim_mol.add_particle(position=center + [0.15, -0.1], color='lightblue', mass=1.0)  # H2

# Add bonds
sim_mol.add_bond(i=0, j=1, k=5000.0, r0=0.15)
sim_mol.add_bond(i=0, j=2, k=5000.0, r0=0.15)
sim_mol.add_angle(i=1, j=0, k=2, k_angle=500.0, theta0=104.5 * np.pi / 180)

# Add more molecules
for i in range(3):
    for j in range(3):
        if i == 0 and j == 0:
            continue
        pos = center + np.array([i * 0.8 - 0.8, j * 0.8 - 0.8])
        o_idx = sim_mol.add_particle(position=pos, color='red', mass=16.0)
        h1_idx = sim_mol.add_particle(position=pos + [0.15, 0.1], color='lightblue', mass=1.0)
        h2_idx = sim_mol.add_particle(position=pos + [0.15, -0.1], color='lightblue', mass=1.0)

        sim_mol.add_bond(i=o_idx, j=h1_idx, k=5000.0, r0=0.15)
        sim_mol.add_bond(i=o_idx, j=h2_idx, k=5000.0, r0=0.15)

sim_mol.initialize_velocities(300.0)

# Equilibrate
for _ in range(500):
    sim_mol.integrate_step()

# Create figure
fig2 = plt.figure(figsize=(16, 10))
gs = GridSpec(2, 3, figure=fig2, wspace=0.3, hspace=0.35)

style_demos = [
    (RenderStyle.SPHERES, ColorScheme.ATOM_TYPE, "Standard Representation"),
    (RenderStyle.BALLS_AND_STICKS, ColorScheme.ATOM_TYPE, "Ball-and-Stick Model"),
    (RenderStyle.SPACE_FILLING, ColorScheme.ATOM_TYPE, "Space-Filling Model"),
    (RenderStyle.HALOS, ColorScheme.ATOM_TYPE, "Glow Effect"),
    (RenderStyle.VELOCITY_ARROWS, ColorScheme.ATOM_TYPE, "With Velocity Vectors"),
    (RenderStyle.WIREFRAME, ColorScheme.ATOM_TYPE, "Wireframe Style"),
]

for idx, (style, color_scheme, title) in enumerate(style_demos):
    row, col = idx // 3, idx % 3
    ax = fig2.add_subplot(gs[row, col])

    renderer = ParticleRenderer(
        style=style,
        color_scheme=color_scheme,
        size_scale=200.0 if style == RenderStyle.SPACE_FILLING else 120.0,
        alpha=0.85
    )

    ax.set_xlim(0, 3.0)
    ax.set_ylim(0, 3.0)
    ax.set_aspect('equal')
    ax.set_title(title, fontweight='bold', pad=10)
    ax.set_xlabel('x (nm)', fontsize=10)
    ax.set_ylabel('y (nm)', fontsize=10)

    # Render bonds first
    if style in [RenderStyle.BALLS_AND_STICKS, RenderStyle.SPACE_FILLING]:
        renderer.render_bonds(ax, sim_mol, bond_thickness=3.0)

    # Render particles
    renderer.render_particles_2d(ax, sim_mol)

    # Panel label
    ax.text(0.02, 0.98, f'({chr(97+idx)})',
           transform=ax.transAxes,
           verticalalignment='top',
           fontsize=14, fontweight='bold',
           bbox=dict(boxstyle='circle', facecolor='lightgreen',
                    edgecolor='black', alpha=0.8))

fig2.suptitle('Visual Representation Styles for Molecular Systems',
             fontsize=15, fontweight='bold', y=0.995)

output2 = os.path.join(os.path.dirname(__file__), '../assets/pub_fig2_styles.png')
fig2.savefig(output2, dpi=300, bbox_inches='tight')
print(f"Saved Figure 2 to: {output2}")

# ===== FIGURE 3: Multi-panel Analysis =====
print("\n" + "=" * 70)
print("Figure 3: Comprehensive System Analysis")
print("=" * 70)

fig3 = plt.figure(figsize=(14, 10))
gs = GridSpec(2, 2, figure=fig3, wspace=0.3, hspace=0.3)

# Use the warm system and run longer simulation
sim_analysis = sim_warm

# Track properties
for _ in range(500):
    sim_analysis.integrate_step()

# Panel (a): Spatial distribution with velocity colors
ax1 = fig3.add_subplot(gs[0, 0])
renderer1 = ParticleRenderer(RenderStyle.HALOS, ColorScheme.VELOCITY,
                             size_scale=180, colormap='plasma')
ax1.set_xlim(0, 3.0)
ax1.set_ylim(0, 3.0)
ax1.set_aspect('equal')
ax1.set_title('(a) Velocity Distribution', fontweight='bold', pad=10)
ax1.set_xlabel('x (nm)')
ax1.set_ylabel('y (nm)')
rect = Rectangle((0, 0), 3.0, 3.0, linewidth=2, edgecolor='black', facecolor='none')
ax1.add_patch(rect)
renderer1.render_particles_2d(ax1, sim_analysis)
create_colorbar(ax1, renderer1, sim_analysis)

# Panel (b): Temperature distribution
ax2 = fig3.add_subplot(gs[0, 1])
renderer2 = ParticleRenderer(RenderStyle.SPHERES, ColorScheme.TEMPERATURE,
                             size_scale=180, colormap='coolwarm')
ax2.set_xlim(0, 3.0)
ax2.set_ylim(0, 3.0)
ax2.set_aspect('equal')
ax2.set_title('(b) Local Temperature', fontweight='bold', pad=10)
ax2.set_xlabel('x (nm)')
ax2.set_ylabel('y (nm)')
rect = Rectangle((0, 0), 3.0, 3.0, linewidth=2, edgecolor='black', facecolor='none')
ax2.add_patch(rect)
renderer2.render_particles_2d(ax2, sim_analysis)
create_colorbar(ax2, renderer2, sim_analysis)

# Panel (c): Energy evolution
ax3 = fig3.add_subplot(gs[1, 0])
times = np.arange(len(sim_analysis.energies['total'])) * sim_analysis.dt
ax3.plot(times, sim_analysis.energies['kinetic'], 'b-', label='Kinetic', linewidth=2)
ax3.plot(times, sim_analysis.energies['potential'], 'r-', label='Potential', linewidth=2)
ax3.plot(times, sim_analysis.energies['total'], 'k-', label='Total', linewidth=2.5)
ax3.set_xlabel('Time (ps)', fontweight='bold')
ax3.set_ylabel('Energy (kJ/mol)', fontweight='bold')
ax3.set_title('(c) Energy Conservation', fontweight='bold', pad=10)
ax3.legend(frameon=True, fancybox=True, shadow=True)
ax3.grid(True, alpha=0.3, linestyle='--')

# Panel (d): Velocity histogram
ax4 = fig3.add_subplot(gs[1, 1])
velocities = sim_analysis.get_velocities()
speeds = np.linalg.norm(velocities, axis=1)
ax4.hist(speeds, bins=20, alpha=0.7, color='steelblue', edgecolor='black', linewidth=1.5)
ax4.axvline(speeds.mean(), color='red', linestyle='--', linewidth=2.5,
           label=f'Mean: {speeds.mean():.3f} nm/ps')
ax4.set_xlabel('Speed (nm/ps)', fontweight='bold')
ax4.set_ylabel('Frequency', fontweight='bold')
ax4.set_title('(d) Velocity Distribution', fontweight='bold', pad=10)
ax4.legend(frameon=True, fancybox=True, shadow=True)
ax4.grid(True, alpha=0.3, linestyle='--', axis='y')

fig3.suptitle(f'Comprehensive Analysis: T = {sim_analysis.temperature():.0f} K',
             fontsize=15, fontweight='bold', y=0.995)

output3 = os.path.join(os.path.dirname(__file__), '../assets/pub_fig3_analysis.png')
fig3.savefig(output3, dpi=300, bbox_inches='tight')
print(f"Saved Figure 3 to: {output3}")

print("\n" + "=" * 70)
print("Publication Figures Complete!")
print("=" * 70)

print(f"""
Created 3 high-quality publication figures:

1. Temperature Effects ({output1})
   - Shows velocity distributions at different temperatures
   - 300 DPI, color-coded by speed
   - Panel labels (a), (b), (c)

2. Rendering Styles ({output2})
   - Demonstrates 6 different visual representations
   - Molecular system with bonds
   - 300 DPI, professional styling

3. Comprehensive Analysis ({output3})
   - Multi-panel analysis figure
   - Spatial + temporal data
   - Energy conservation + statistics
   - 300 DPI, publication-ready

All figures suitable for:
- Research papers
- Presentations
- Educational materials
- Posters and slides

Tips for publication:
- Use 300 DPI for print, 150 DPI for web
- Vector formats (PDF/SVG) for scalability
- Consistent color schemes across figures
- Clear labels and legends
- High contrast for readability
""")

plt.show()
