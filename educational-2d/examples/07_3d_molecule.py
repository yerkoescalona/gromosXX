"""
Example 7: 3D Molecule (Water-like Structure)
==============================================

Demonstrates:
- Bonded system in 3D
- Angular constraints creating molecular shape
- Molecular vibrations and rotations
- 3D visualization of molecular structure

This creates a simple water-like molecule (3 atoms in bent geometry)
and watches it vibrate and rotate in 3D space.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from core.simulation3d import Simulation3D
from core.renderer3d import Renderer3D

# Simulation parameters
BOX_SIZE = (3.0, 3.0, 3.0)
TEMPERATURE = 300.0  # K
BOND_LENGTH = 0.15  # nm
BOND_K = 5000.0  # stiff bonds
ANGLE_K = 500.0  # stiff angle (to maintain bent shape)
ANGLE_0 = 104.5 * np.pi / 180  # Water-like angle (~104.5 degrees)

print("=" * 70)
print("Educational 3D Simulation: Water-like Molecule")
print("=" * 70)
print(f"\nMolecule Parameters:")
print(f"  Structure: H-O-H (bent)")
print(f"  Bond length: {BOND_LENGTH} nm")
print(f"  H-O-H angle: {ANGLE_0 * 180 / np.pi:.1f}°")
print(f"  Temperature: {TEMPERATURE} K")
print(f"  No periodic boundaries (single molecule in vacuum)")
print()

# Create simulation
sim = Simulation3D(
    box_size=BOX_SIZE,
    periodic=False,  # No PBC for isolated molecule
    integrator='velocity-verlet',
    dt=0.0005  # Smaller timestep for stiff bonds
)

# Build water-like molecule
print("Building molecule...")

# Center of box
center = np.array(BOX_SIZE) / 2

# Oxygen atom (central, red)
o_idx = sim.add_particle(
    position=center,
    mass=16.0,  # Oxygen mass
    color='red'
)

# Hydrogen 1 (blue)
angle = ANGLE_0 / 2
h1_pos = center + np.array([BOND_LENGTH * np.cos(angle), BOND_LENGTH * np.sin(angle), 0])
h1_idx = sim.add_particle(
    position=h1_pos,
    mass=1.0,  # Hydrogen mass
    color='lightblue'
)

# Hydrogen 2 (blue)
h2_pos = center + np.array([BOND_LENGTH * np.cos(angle), -BOND_LENGTH * np.sin(angle), 0])
h2_idx = sim.add_particle(
    position=h2_pos,
    mass=1.0,
    color='lightblue'
)

# Add bonds: H1-O and H2-O
sim.add_bond(i=h1_idx, j=o_idx, k=BOND_K, r0=BOND_LENGTH)
sim.add_bond(i=h2_idx, j=o_idx, k=BOND_K, r0=BOND_LENGTH)

# Add angle: H1-O-H2
sim.add_angle(i=h1_idx, j=o_idx, k=h2_idx, k_angle=ANGLE_K, theta0=ANGLE_0)

print(f"Created molecule:")
print(f"  Atoms: {len(sim.particles)} (1 O + 2 H)")
print(f"  Bonds: {len(sim.bonds)}")
print(f"  Angles: {len(sim.angles)}")

# Initialize velocities
print(f"\nInitializing velocities at {TEMPERATURE} K...")
sim.initialize_velocities(TEMPERATURE)

# Light temperature coupling
sim.set_temperature(TEMPERATURE, tau=0.5)

print(f"Initial temperature: {sim.temperature():.1f} K")

# Create renderer
renderer = Renderer3D(
    box_size=BOX_SIZE,
    figsize=(14, 6),
    azimuth=45,
    elevation=20
)

# Analyze molecular geometry
def analyze_molecule(sim):
    """Calculate molecular properties."""
    o_pos = sim.particles[0].position
    h1_pos = sim.particles[1].position
    h2_pos = sim.particles[2].position

    # Bond lengths
    r_oh1 = np.linalg.norm(h1_pos - o_pos)
    r_oh2 = np.linalg.norm(h2_pos - o_pos)

    # H-O-H angle
    v1 = h1_pos - o_pos
    v2 = h2_pos - o_pos
    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.arccos(np.clip(cos_angle, -1, 1))

    # Dipole moment direction (simplified)
    h_avg = (h1_pos + h2_pos) / 2
    dipole_dir = o_pos - h_avg

    return r_oh1, r_oh2, angle, dipole_dir

# Initial geometry
r_oh1, r_oh2, angle, _ = analyze_molecule(sim)
print(f"\nInitial Geometry:")
print(f"  O-H1 bond: {r_oh1:.4f} nm")
print(f"  O-H2 bond: {r_oh2:.4f} nm")
print(f"  H-O-H angle: {angle * 180 / np.pi:.2f}°")

# Run simulation
print("\n" + "=" * 70)
print("Simulation: Molecular Vibrations and Rotations")
print("=" * 70)
print("\nObserve:")
print("  - Bond stretching vibrations")
print("  - Angle bending vibrations")
print("  - Overall molecular rotation")
print("  - Molecule tumbling in space")
print()

# Track molecular properties
bond_lengths = []
angles_tracked = []

for step in range(3000):
    sim.integrate_step()

    # Track geometry
    r_oh1, r_oh2, angle, _ = analyze_molecule(sim)
    bond_lengths.append((r_oh1 + r_oh2) / 2)
    angles_tracked.append(angle * 180 / np.pi)

    # Render every 10 steps with rotating view
    if step % 10 == 0:
        azimuth = 45 + (720 * step / 3000)  # Two full rotations
        renderer.render(sim, show_bonds=True, azimuth=azimuth, elevation=20)

    if step % 600 == 0:
        ke, pe, total = sim.total_energy()
        print(f"Step {step:4d}: T = {sim.temperature():5.1f} K, "
              f"Bond = {(r_oh1+r_oh2)/2:.4f} nm, "
              f"Angle = {angle*180/np.pi:6.2f}°, "
              f"E = {total:7.2f} kJ/mol")

print("\n" + "=" * 70)
print("Simulation Complete!")
print("=" * 70)

# Final analysis
r_oh1, r_oh2, angle, _ = analyze_molecule(sim)
ke, pe, total = sim.total_energy()

print(f"\nFinal Molecular Geometry:")
print(f"  O-H1 bond: {r_oh1:.4f} nm (equilibrium: {BOND_LENGTH:.4f} nm)")
print(f"  O-H2 bond: {r_oh2:.4f} nm (equilibrium: {BOND_LENGTH:.4f} nm)")
print(f"  H-O-H angle: {angle * 180 / np.pi:.2f}° (equilibrium: {ANGLE_0 * 180 / np.pi:.1f}°)")

print(f"\nEnergy Analysis:")
print(f"  Total energy: {total:.2f} kJ/mol")
print(f"  Bond energy: {sim.energies['bond'][-1]:.2f} kJ/mol")
print(f"  Angle energy: {sim.energies['angle'][-1]:.2f} kJ/mol")
print(f"  Temperature: {sim.temperature():.1f} K")

# Vibrational analysis
bond_mean = np.mean(bond_lengths)
bond_std = np.std(bond_lengths)
angle_mean = np.mean(angles_tracked)
angle_std = np.std(angles_tracked)

print(f"\nVibrational Statistics:")
print(f"  Bond length: {bond_mean:.4f} ± {bond_std:.4f} nm")
print(f"  H-O-H angle: {angle_mean:.2f} ± {angle_std:.2f}°")

# Plot vibrational motion
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(12, 8))

times = np.arange(len(bond_lengths)) * sim.dt

# Bond length vibrations
axes[0, 0].plot(times, bond_lengths, 'b-', linewidth=0.5, alpha=0.7)
axes[0, 0].axhline(BOND_LENGTH, color='r', linestyle='--', label=f'Equilibrium: {BOND_LENGTH:.3f} nm')
axes[0, 0].set_xlabel('Time (ps)')
axes[0, 0].set_ylabel('Average O-H Bond Length (nm)')
axes[0, 0].set_title('Bond Stretching Vibration')
axes[0, 0].legend()
axes[0, 0].grid(True, alpha=0.3)

# Angle vibrations
axes[0, 1].plot(times, angles_tracked, 'g-', linewidth=0.5, alpha=0.7)
axes[0, 1].axhline(ANGLE_0 * 180 / np.pi, color='r', linestyle='--',
                   label=f'Equilibrium: {ANGLE_0*180/np.pi:.1f}°')
axes[0, 1].set_xlabel('Time (ps)')
axes[0, 1].set_ylabel('H-O-H Angle (degrees)')
axes[0, 1].set_title('Angle Bending Vibration')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.3)

# Bond length distribution
axes[1, 0].hist(bond_lengths, bins=50, alpha=0.7, edgecolor='black')
axes[1, 0].axvline(BOND_LENGTH, color='r', linestyle='--', linewidth=2,
                   label=f'Equilibrium: {BOND_LENGTH:.3f} nm')
axes[1, 0].axvline(bond_mean, color='b', linestyle='--', linewidth=2,
                   label=f'Mean: {bond_mean:.3f} nm')
axes[1, 0].set_xlabel('Bond Length (nm)')
axes[1, 0].set_ylabel('Frequency')
axes[1, 0].set_title('Bond Length Distribution')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3)

# Angle distribution
axes[1, 1].hist(angles_tracked, bins=50, alpha=0.7, edgecolor='black', color='green')
axes[1, 1].axvline(ANGLE_0 * 180 / np.pi, color='r', linestyle='--', linewidth=2,
                   label=f'Equilibrium: {ANGLE_0*180/np.pi:.1f}°')
axes[1, 1].axvline(angle_mean, color='b', linestyle='--', linewidth=2,
                   label=f'Mean: {angle_mean:.1f}°')
axes[1, 1].set_xlabel('H-O-H Angle (degrees)')
axes[1, 1].set_ylabel('Frequency')
axes[1, 1].set_title('Angle Distribution')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()

vibration_file = os.path.join(os.path.dirname(__file__), '../assets/molecule_3d_vibrations.png')
fig.savefig(vibration_file, dpi=150, bbox_inches='tight')
print(f"\nSaved vibrational analysis to: {vibration_file}")

# Save final molecular structure
output_file = os.path.join(os.path.dirname(__file__), '../assets/molecule_3d.png')
renderer.save_frame(output_file)
print(f"Saved molecular structure to: {output_file}")

# Try to create rotating movie
try:
    from core.movie_utils import make_rotating_3d_movie

    output_movie = os.path.join(os.path.dirname(__file__), '../assets/molecule_3d_rotation.mp4')
    print(f"\nCreating rotating molecule movie...")

    make_rotating_3d_movie(
        sim, renderer,
        n_frames=360,
        output_file=output_movie,
        elevation=20,
        fps=30
    )
    print(f"Saved rotating movie to: {output_movie}")

except Exception as e:
    print(f"\nCould not create movie: {e}")

print("\n" + "=" * 70)
print("Educational Notes:")
print("=" * 70)
print("""
Molecular Vibrations in 3D:

1. Normal Modes:
   This water-like molecule has normal vibrational modes:
   - Symmetric stretch: Both O-H bonds stretch together
   - Asymmetric stretch: One O-H stretches, other compresses
   - Bending: H-O-H angle changes

2. Thermal Motion:
   At finite temperature, all modes are excited
   Classical limit: Each mode gets k_B*T/2 energy on average

3. Rotation vs Vibration:
   - Vibrations: Internal coordinates change (bonds, angles)
   - Rotations: Whole molecule tumbles in space
   - Both occur simultaneously at T > 0

4. Comparison with Real Water:
   - Real H2O has quantum vibrational levels
   - Real bonds are anharmonic (not purely harmonic)
   - Real H2O has electrostatics (we simplified to LJ + bonds)

5. This simulation shows:
   - Classical molecular vibrations
   - Thermal averaging
   - Coupling between rotation and vibration
   - Energy exchange between modes
""")

renderer.show()
plt.show()
