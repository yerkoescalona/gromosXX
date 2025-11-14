"""
Advanced Rendering Styles for Molecular Visualization
======================================================

Multiple visual representations for particles/molecules:
- Size scaling (by mass, charge, temperature, etc.)
- Color schemes (velocity, energy, atom type, etc.)
- Visual styles (spheres, points, halos, etc.)
- Transparency and opacity effects
- Ball-and-stick vs space-filling models
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from typing import List, Optional, Callable, Tuple
from enum import Enum


class RenderStyle(Enum):
    """Available rendering styles."""
    SPHERES = "spheres"              # Standard spheres
    POINTS = "points"                # Simple points
    HALOS = "halos"                  # Spheres with glow
    VDW_RADII = "vdw"               # Van der Waals radii
    BALLS_AND_STICKS = "ball_stick" # Ball-and-stick model
    SPACE_FILLING = "space_fill"    # Space-filling (CPK)
    WIREFRAME = "wireframe"          # Wireframe spheres
    VELOCITY_ARROWS = "vel_arrows"  # Velocity vectors


class ColorScheme(Enum):
    """Available color schemes."""
    ATOM_TYPE = "atom_type"          # By atom type (fixed colors)
    VELOCITY = "velocity"            # By velocity magnitude
    KINETIC_ENERGY = "kinetic"       # By kinetic energy
    TEMPERATURE = "temperature"      # Local temperature
    FORCE = "force"                  # By force magnitude
    CUSTOM = "custom"                # User-provided colors
    CPK = "cpk"                      # Standard CPK colors
    RAINBOW = "rainbow"              # Rainbow by index


class ParticleRenderer:
    """
    Advanced particle rendering with multiple visual styles.

    Supports:
    - Different rendering styles (spheres, halos, etc.)
    - Multiple color schemes
    - Size scaling
    - Transparency effects
    - Custom colormaps
    """

    def __init__(self,
                 style: RenderStyle = RenderStyle.SPHERES,
                 color_scheme: ColorScheme = ColorScheme.ATOM_TYPE,
                 size_scale: float = 100.0,
                 alpha: float = 0.8,
                 colormap: str = 'viridis'):
        """
        Initialize particle renderer.

        Parameters
        ----------
        style : RenderStyle
            Visual representation style
        color_scheme : ColorScheme
            Color coding scheme
        size_scale : float
            Base size for particles
        alpha : float
            Transparency (0=transparent, 1=opaque)
        colormap : str
            Matplotlib colormap name for continuous schemes
        """
        self.style = style
        self.color_scheme = color_scheme
        self.size_scale = size_scale
        self.alpha = alpha
        self.colormap = colormap

        # CPK color scheme (standard molecular colors)
        self.cpk_colors = {
            'H': '#FFFFFF',   # Hydrogen - white
            'C': '#909090',   # Carbon - gray
            'N': '#3050F8',   # Nitrogen - blue
            'O': '#FF0D0D',   # Oxygen - red
            'F': '#90E050',   # Fluorine - green
            'P': '#FF8000',   # Phosphorus - orange
            'S': '#FFFF30',   # Sulfur - yellow
            'Cl': '#1FF01F',  # Chlorine - green
            'Default': '#FF1493'  # Magenta for unknown
        }

    def get_colors(self, sim, particles_subset: Optional[List[int]] = None):
        """
        Get colors for particles based on color scheme.

        Parameters
        ----------
        sim : Simulation2D or Simulation3D
            Simulation object
        particles_subset : List[int], optional
            Indices of particles to render (None = all)

        Returns
        -------
        colors : array-like
            Colors for each particle
        """
        if particles_subset is None:
            particles = sim.particles
        else:
            particles = [sim.particles[i] for i in particles_subset]

        if self.color_scheme == ColorScheme.ATOM_TYPE:
            # Use particle's color attribute
            return [p.color for p in particles]

        elif self.color_scheme == ColorScheme.VELOCITY:
            # Color by velocity magnitude
            velocities = np.array([p.velocity for p in particles])
            speeds = np.linalg.norm(velocities, axis=1)
            return self._map_to_colormap(speeds)

        elif self.color_scheme == ColorScheme.KINETIC_ENERGY:
            # Color by kinetic energy
            energies = np.array([p.kinetic_energy() for p in particles])
            return self._map_to_colormap(energies)

        elif self.color_scheme == ColorScheme.FORCE:
            # Color by force magnitude
            forces = np.array([p.force for p in particles])
            force_mags = np.linalg.norm(forces, axis=1)
            return self._map_to_colormap(force_mags)

        elif self.color_scheme == ColorScheme.TEMPERATURE:
            # Color by local kinetic temperature
            # T_local = 2*KE / (k_B * DOF)
            k_B = 0.00831446261815324
            dof = len(particles[0].velocity)
            temps = np.array([2 * p.kinetic_energy() / (k_B * dof)
                            for p in particles])
            return self._map_to_colormap(temps)

        elif self.color_scheme == ColorScheme.CPK:
            # Standard CPK colors (would need atom type info)
            return [self.cpk_colors.get(p.color, self.cpk_colors['Default'])
                   for p in particles]

        elif self.color_scheme == ColorScheme.RAINBOW:
            # Rainbow by particle index
            indices = np.arange(len(particles))
            return self._map_to_colormap(indices)

        else:  # CUSTOM
            return [p.color for p in particles]

    def get_sizes(self, sim, particles_subset: Optional[List[int]] = None):
        """
        Get sizes for particles based on style and properties.

        Parameters
        ----------
        sim : Simulation2D or Simulation3D
            Simulation object
        particles_subset : List[int], optional
            Indices of particles to render

        Returns
        -------
        sizes : array
            Sizes for each particle
        """
        if particles_subset is None:
            particles = sim.particles
        else:
            particles = [sim.particles[i] for i in particles_subset]

        if self.style == RenderStyle.VDW_RADII:
            # Size by Van der Waals radius (from sigma parameter)
            sizes = []
            for p in particles:
                if hasattr(sim, 'lj_params') and p.atom_type in sim.lj_params:
                    sigma = sim.lj_params[p.atom_type].get('sigma', 0.3)
                    # Convert sigma to display size (arbitrary scaling)
                    sizes.append(self.size_scale * (sigma / 0.3) ** 2)
                else:
                    sizes.append(self.size_scale)
            return np.array(sizes)

        elif self.style == RenderStyle.SPACE_FILLING:
            # Larger spheres that touch
            return np.ones(len(particles)) * self.size_scale * 2.0

        elif self.style == RenderStyle.BALLS_AND_STICKS:
            # Smaller spheres for ball-and-stick
            return np.ones(len(particles)) * self.size_scale * 0.5

        elif self.style == RenderStyle.POINTS:
            # Very small points
            return np.ones(len(particles)) * self.size_scale * 0.1

        else:  # Default spheres
            return np.ones(len(particles)) * self.size_scale

    def render_particles_2d(self, ax, sim, particles_subset: Optional[List[int]] = None):
        """
        Render particles in 2D with selected style.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to render on
        sim : Simulation2D
            Simulation object
        particles_subset : List[int], optional
            Particle indices to render
        """
        positions = sim.get_positions()

        if particles_subset is not None:
            positions = positions[particles_subset]

        colors = self.get_colors(sim, particles_subset)
        sizes = self.get_sizes(sim, particles_subset)

        if self.style == RenderStyle.HALOS:
            # Draw halos (larger transparent outer circle)
            ax.scatter(positions[:, 0], positions[:, 1],
                      c=colors, s=sizes * 3, alpha=self.alpha * 0.3,
                      edgecolors='none')
            # Draw main spheres
            ax.scatter(positions[:, 0], positions[:, 1],
                      c=colors, s=sizes, alpha=self.alpha,
                      edgecolors='black', linewidth=0.5)

        elif self.style == RenderStyle.WIREFRAME:
            # Wireframe circles
            ax.scatter(positions[:, 0], positions[:, 1],
                      c='none', s=sizes, alpha=1.0,
                      edgecolors=colors, linewidth=2)

        elif self.style == RenderStyle.VELOCITY_ARROWS:
            # Draw particles
            ax.scatter(positions[:, 0], positions[:, 1],
                      c=colors, s=sizes, alpha=self.alpha,
                      edgecolors='black', linewidth=0.5)
            # Draw velocity arrows
            velocities = sim.get_velocities()
            if particles_subset is not None:
                velocities = velocities[particles_subset]
            scale = 0.05
            for pos, vel in zip(positions, velocities):
                if np.linalg.norm(vel) > 0.01:
                    ax.arrow(pos[0], pos[1], vel[0]*scale, vel[1]*scale,
                           head_width=0.03, head_length=0.03,
                           fc='red', ec='red', alpha=0.6, width=0.01)

        else:  # Standard spheres
            ax.scatter(positions[:, 0], positions[:, 1],
                      c=colors, s=sizes, alpha=self.alpha,
                      edgecolors='black', linewidth=1)

    def render_particles_3d(self, ax, sim, particles_subset: Optional[List[int]] = None):
        """
        Render particles in 3D with selected style.

        Parameters
        ----------
        ax : mpl_toolkits.mplot3d.axes3d.Axes3D
            3D axes to render on
        sim : Simulation3D
            Simulation object
        particles_subset : List[int], optional
            Particle indices to render
        """
        positions = sim.get_positions()

        if particles_subset is not None:
            positions = positions[particles_subset]

        colors = self.get_colors(sim, particles_subset)
        sizes = self.get_sizes(sim, particles_subset)

        if self.style == RenderStyle.HALOS:
            # Outer glow
            ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
                      c=colors, s=sizes * 3, alpha=self.alpha * 0.3,
                      edgecolors='none')
            # Main spheres
            ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
                      c=colors, s=sizes, alpha=self.alpha,
                      edgecolors='black', linewidth=0.5)

        elif self.style == RenderStyle.WIREFRAME:
            # Wireframe spheres
            ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
                      c='none', s=sizes, alpha=1.0,
                      edgecolors=colors, linewidth=2)

        else:  # Standard
            ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
                      c=colors, s=sizes, alpha=self.alpha,
                      edgecolors='black', linewidth=1)

    def render_bonds(self, ax, sim, bond_thickness: float = 2.0,
                    bond_color: str = 'gray', is_3d: bool = False):
        """
        Render bonds between particles.

        Parameters
        ----------
        ax : matplotlib.axes.Axes or Axes3D
            Axes to render on
        sim : Simulation2D or Simulation3D
            Simulation object
        bond_thickness : float
            Line thickness for bonds
        bond_color : str
            Bond color
        is_3d : bool
            Whether rendering in 3D
        """
        positions = sim.get_positions()

        for bond in sim.bonds:
            pi = positions[bond.i]
            pj = positions[bond.j]

            # Handle periodic boundaries
            dr = pj - pi
            if sim.periodic:
                dr = dr - sim.box_size * np.round(dr / sim.box_size)
                pj_wrapped = pi + dr

                if is_3d:
                    ax.plot3D([pi[0], pj_wrapped[0]],
                             [pi[1], pj_wrapped[1]],
                             [pi[2], pj_wrapped[2]],
                             color=bond_color, linewidth=bond_thickness,
                             alpha=0.6)
                else:
                    ax.plot([pi[0], pj_wrapped[0]],
                           [pi[1], pj_wrapped[1]],
                           color=bond_color, linewidth=bond_thickness,
                           alpha=0.6)
            else:
                if is_3d:
                    ax.plot3D([pi[0], pj[0]], [pi[1], pj[1]], [pi[2], pj[2]],
                             color=bond_color, linewidth=bond_thickness,
                             alpha=0.6)
                else:
                    ax.plot([pi[0], pj[0]], [pi[1], pj[1]],
                           color=bond_color, linewidth=bond_thickness,
                           alpha=0.6)

    def _map_to_colormap(self, values: np.ndarray):
        """Map values to colormap colors."""
        if len(values) == 0:
            return []

        # Normalize values to [0, 1]
        vmin, vmax = values.min(), values.max()
        if vmax - vmin < 1e-10:
            normalized = np.ones_like(values) * 0.5
        else:
            normalized = (values - vmin) / (vmax - vmin)

        # Get colormap
        cmap = plt.get_cmap(self.colormap)

        # Map to colors
        return [cmap(val) for val in normalized]


def create_colorbar(ax, renderer: ParticleRenderer, sim,
                   label: str = "Property"):
    """
    Add a colorbar for continuous color schemes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to add colorbar to
    renderer : ParticleRenderer
        Renderer with color scheme
    sim : Simulation2D or Simulation3D
        Simulation object
    label : str
        Colorbar label
    """
    if renderer.color_scheme in [ColorScheme.VELOCITY, ColorScheme.KINETIC_ENERGY,
                                  ColorScheme.FORCE, ColorScheme.TEMPERATURE]:
        # Get data range
        particles = sim.particles

        if renderer.color_scheme == ColorScheme.VELOCITY:
            velocities = np.array([p.velocity for p in particles])
            values = np.linalg.norm(velocities, axis=1)
            label = "Velocity (nm/ps)"
        elif renderer.color_scheme == ColorScheme.KINETIC_ENERGY:
            values = np.array([p.kinetic_energy() for p in particles])
            label = "Kinetic Energy (kJ/mol)"
        elif renderer.color_scheme == ColorScheme.FORCE:
            forces = np.array([p.force for p in particles])
            values = np.linalg.norm(forces, axis=1)
            label = "Force Magnitude (kJ/mol/nm)"
        else:  # TEMPERATURE
            k_B = 0.00831446261815324
            dof = len(particles[0].velocity)
            values = np.array([2 * p.kinetic_energy() / (k_B * dof)
                             for p in particles])
            label = "Temperature (K)"

        # Create colorbar
        sm = plt.cm.ScalarMappable(cmap=renderer.colormap,
                                   norm=plt.Normalize(vmin=values.min(),
                                                     vmax=values.max()))
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label=label, fraction=0.046, pad=0.04)


# Example usage and demonstrations
if __name__ == "__main__":
    print("Advanced Particle Rendering Styles")
    print("=" * 60)
    print("\nAvailable styles:")
    for style in RenderStyle:
        print(f"  - {style.value}")

    print("\nAvailable color schemes:")
    for scheme in ColorScheme:
        print(f"  - {scheme.value}")

    print("\nUsage example:")
    print("""
from core.render_styles import ParticleRenderer, RenderStyle, ColorScheme

# Create renderer with specific style
renderer = ParticleRenderer(
    style=RenderStyle.HALOS,
    color_scheme=ColorScheme.VELOCITY,
    size_scale=150.0,
    alpha=0.8,
    colormap='plasma'
)

# Render particles
fig, ax = plt.subplots()
renderer.render_particles_2d(ax, sim)
renderer.render_bonds(ax, sim)
    """)
