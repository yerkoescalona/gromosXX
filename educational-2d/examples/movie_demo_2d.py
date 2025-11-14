"""
Movie Demo: Creating Animations from 2D Simulations
====================================================

Demonstrates how to easily create movies from 2D MD simulations.

Supports:
- MP4 videos
- Animated GIFs
- Frame sequences

Requirements:
- ffmpeg (for MP4): conda install -c conda-forge ffmpeg
- Pillow (for GIF): pip install Pillow
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from core import Simulation2D, Renderer2D
from core.movie_utils import (
    make_quick_movie,
    make_high_quality_movie,
    make_gif_animation,
    record_simulation_movie
)

print("=" * 70)
print("Educational 2D Simulation: Movie Creation Demo")
print("=" * 70)

# Create a simple gas simulation
N_PARTICLES = 30
BOX_SIZE = (4.0, 4.0)
TEMPERATURE = 300.0

print(f"\nCreating simulation...")
print(f"  Particles: {N_PARTICLES}")
print(f"  Box: {BOX_SIZE[0]} x {BOX_SIZE[1]} nm")
print(f"  Temperature: {TEMPERATURE} K")

sim = Simulation2D(
    box_size=BOX_SIZE,
    periodic=True,
    integrator='velocity-verlet',
    dt=0.002
)

# Add particles in grid
grid_size = int(np.ceil(np.sqrt(N_PARTICLES)))
spacing = min(BOX_SIZE) / grid_size

for i in range(N_PARTICLES):
    ix = i % grid_size
    iy = i // grid_size
    x = (ix + 0.5) * spacing + np.random.uniform(-0.05, 0.05)
    y = (iy + 0.5) * spacing + np.random.uniform(-0.05, 0.05)

    color = ['blue', 'red', 'green', 'orange'][i % 4]
    sim.add_particle(position=[x, y], color=color)

# Initialize velocities
sim.initialize_velocities(TEMPERATURE)
sim.set_temperature(TEMPERATURE, tau=0.1)

# Create renderer
renderer = Renderer2D(
    box_size=BOX_SIZE,
    figsize=(14, 6),
    show_trails=True,
    trail_length=30
)

# Equilibrate briefly
print("\nEquilibrating (500 steps)...")
for step in range(500):
    sim.integrate_step()
    if step % 100 == 0:
        print(f"  Step {step}/500")

print("\nEquilibration complete!")

# --- Movie Option 1: Quick Movie (MP4) ---
print("\n" + "=" * 70)
print("Option 1: Quick Movie (300 frames, 30 FPS)")
print("=" * 70)

output_dir = os.path.join(os.path.dirname(__file__), '../assets')
os.makedirs(output_dir, exist_ok=True)

try:
    make_quick_movie(
        sim, renderer,
        n_steps=1000,
        filename=os.path.join(output_dir, 'gas_quick.mp4')
    )
except Exception as e:
    print(f"Could not create MP4: {e}")
    print("Install ffmpeg: conda install -c conda-forge ffmpeg")

# Reset simulation for next demo
sim.step = 0
sim.time = 0.0
sim.energies = {key: [] for key in sim.energies.keys()}

# --- Movie Option 2: Animated GIF ---
print("\n" + "=" * 70)
print("Option 2: Animated GIF (smaller file, web-friendly)")
print("=" * 70)

try:
    make_gif_animation(
        sim, renderer,
        n_steps=500,
        filename=os.path.join(output_dir, 'gas_animation.gif')
    )
except Exception as e:
    print(f"Could not create GIF: {e}")
    print("Install Pillow: pip install Pillow")

# Reset again
sim.step = 0
sim.time = 0.0
sim.energies = {key: [] for key in sim.energies.keys()}

# --- Movie Option 3: High Quality (every frame) ---
print("\n" + "=" * 70)
print("Option 3: High Quality Movie (60 FPS, high DPI)")
print("=" * 70)
print("This will be slower but higher quality...")

try:
    # Shorter for demo
    make_high_quality_movie(
        sim, renderer,
        n_steps=300,
        filename=os.path.join(output_dir, 'gas_highquality.mp4')
    )
except Exception as e:
    print(f"Could not create high-quality MP4: {e}")

# --- Movie Option 4: Custom Settings ---
print("\n" + "=" * 70)
print("Option 4: Custom Movie Settings")
print("=" * 70)

sim.step = 0
sim.time = 0.0
sim.energies = {key: [] for key in sim.energies.keys()}

try:
    record_simulation_movie(
        sim, renderer,
        n_steps=800,
        output_file=os.path.join(output_dir, 'gas_custom.mp4'),
        render_interval=2,  # Render every 2 steps
        fps=24,            # 24 FPS (cinematic)
        dpi=120,           # Medium-high resolution
        show_bonds=False,
        progress_interval=100
    )
except Exception as e:
    print(f"Could not create custom movie: {e}")

print("\n" + "=" * 70)
print("Movie Creation Demo Complete!")
print("=" * 70)
print(f"\nCheck the output directory: {output_dir}")
print("\nFiles created:")
print("  - gas_quick.mp4 (quick movie)")
print("  - gas_animation.gif (animated GIF)")
print("  - gas_highquality.mp4 (high quality)")
print("  - gas_custom.mp4 (custom settings)")
print("\nTips:")
print("  - MP4 files: Good for presentations and video players")
print("  - GIF files: Good for websites and documentation")
print("  - Adjust FPS and DPI for quality vs file size tradeoff")
print("  - Use render_interval to skip frames for faster rendering")
