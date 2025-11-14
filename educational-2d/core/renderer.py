"""
2D Visualization and Rendering
===============================

Real-time matplotlib visualization for 2D molecular dynamics simulations.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
from matplotlib.collections import LineCollection
from typing import Optional, List, Tuple


class Renderer2D:
    """
    Real-time visualization for 2D MD simulations.

    Features:
    - Particle rendering with colors
    - Bond visualization
    - Force vectors
    - Energy plots
    - Temperature display
    - Trajectory trails
    """

    def __init__(self,
                 box_size: Tuple[float, float],
                 figsize: Tuple[float, float] = (12, 6),
                 show_forces: bool = False,
                 show_velocities: bool = False,
                 show_trails: bool = False,
                 trail_length: int = 50):
        """
        Initialize 2D renderer.

        Parameters
        ----------
        box_size : Tuple[float, float]
            Simulation box size (width, height) in nm
        figsize : Tuple[float, float]
            Figure size in inches
        show_forces : bool
            Display force vectors on particles
        show_velocities : bool
            Display velocity vectors on particles
        show_trails : bool
            Show particle trajectories
        trail_length : int
            Number of previous positions to show in trail
        """
        self.box_size = np.array(box_size)
        self.figsize = figsize
        self.show_forces = show_forces
        self.show_velocities = show_velocities
        self.show_trails = show_trails
        self.trail_length = trail_length

        # Trajectory storage for trails
        self.trails = []

        # Create figure with subplots
        self.fig, self.axes = plt.subplots(1, 2, figsize=figsize)
        self.ax_sim = self.axes[0]   # Simulation view
        self.ax_energy = self.axes[1]  # Energy plot

        # Setup simulation view
        self.ax_sim.set_xlim(0, box_size[0])
        self.ax_sim.set_ylim(0, box_size[1])
        self.ax_sim.set_aspect('equal')
        self.ax_sim.set_xlabel('x (nm)')
        self.ax_sim.set_ylabel('y (nm)')
        self.ax_sim.set_title('2D Molecular Dynamics')

        # Draw box
        rect = patches.Rectangle((0, 0), box_size[0], box_size[1],
                                 linewidth=2, edgecolor='black',
                                 facecolor='none')
        self.ax_sim.add_patch(rect)

        # Setup energy plot
        self.ax_energy.set_xlabel('Time (ps)')
        self.ax_energy.set_ylabel('Energy (kJ/mol)')
        self.ax_energy.set_title('Energy Evolution')
        self.ax_energy.grid(True, alpha=0.3)

        # Plot objects (will be initialized later)
        self.scatter = None
        self.bond_lines = None
        self.force_quiver = None
        self.velocity_quiver = None
        self.energy_lines = {}
        self.text_info = None

        plt.tight_layout()

    def render_frame(self, sim, show_bonds: bool = True):
        """
        Render a single frame of the simulation.

        Parameters
        ----------
        sim : Simulation2D
            The simulation object to render
        show_bonds : bool
            Whether to display bonds
        """
        positions = sim.get_positions()
        colors = sim.get_colors()

        # Clear previous frame
        self.ax_sim.clear()

        # Redraw box
        self.ax_sim.set_xlim(0, self.box_size[0])
        self.ax_sim.set_ylim(0, self.box_size[1])
        self.ax_sim.set_aspect('equal')
        self.ax_sim.set_xlabel('x (nm)')
        self.ax_sim.set_ylabel('y (nm)')

        rect = patches.Rectangle((0, 0), self.box_size[0], self.box_size[1],
                                 linewidth=2, edgecolor='black',
                                 facecolor='none')
        self.ax_sim.add_patch(rect)

        # Draw trails
        if self.show_trails and len(self.trails) > 0:
            for particle_idx in range(len(positions)):
                if particle_idx < len(self.trails):
                    trail = np.array(self.trails[particle_idx])
                    if len(trail) > 1:
                        self.ax_sim.plot(trail[:, 0], trail[:, 1],
                                       'k-', alpha=0.2, linewidth=0.5)

        # Draw bonds
        if show_bonds and len(sim.bonds) > 0:
            for bond in sim.bonds:
                pi = positions[bond.i]
                pj = positions[bond.j]

                # Handle periodic boundaries for bond visualization
                dr = pj - pi
                if sim.periodic:
                    dr = dr - self.box_size * np.round(dr / self.box_size)
                    pj_wrapped = pi + dr

                    self.ax_sim.plot([pi[0], pj_wrapped[0]],
                                    [pi[1], pj_wrapped[1]],
                                    'gray', linewidth=2, alpha=0.6)
                else:
                    self.ax_sim.plot([pi[0], pj[0]], [pi[1], pj[1]],
                                    'gray', linewidth=2, alpha=0.6)

        # Draw particles
        self.ax_sim.scatter(positions[:, 0], positions[:, 1],
                          c=colors, s=100, alpha=0.8, edgecolors='black',
                          linewidth=1)

        # Draw force vectors
        if self.show_forces:
            forces = sim.get_forces()
            # Scale forces for visualization
            force_scale = 0.05
            for i, (pos, force) in enumerate(zip(positions, forces)):
                if np.linalg.norm(force) > 0.01:
                    self.ax_sim.arrow(pos[0], pos[1],
                                    force[0] * force_scale,
                                    force[1] * force_scale,
                                    head_width=0.05, head_length=0.05,
                                    fc='red', ec='red', alpha=0.6)

        # Draw velocity vectors
        if self.show_velocities:
            velocities = sim.get_velocities()
            vel_scale = 0.1
            for i, (pos, vel) in enumerate(zip(positions, velocities)):
                if np.linalg.norm(vel) > 0.01:
                    self.ax_sim.arrow(pos[0], pos[1],
                                    vel[0] * vel_scale,
                                    vel[1] * vel_scale,
                                    head_width=0.05, head_length=0.05,
                                    fc='blue', ec='blue', alpha=0.6)

        # Display simulation info
        info_text = f"Step: {sim.step}\n"
        info_text += f"Time: {sim.time:.3f} ps\n"
        info_text += f"T: {sim.temperature():.1f} K\n"
        info_text += f"Particles: {len(sim.particles)}"

        self.ax_sim.text(0.02, 0.98, info_text,
                        transform=self.ax_sim.transAxes,
                        verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                        fontsize=9)

        self.ax_sim.set_title('2D Molecular Dynamics')

    def update_energy_plot(self, sim):
        """
        Update the energy evolution plot.

        Parameters
        ----------
        sim : Simulation2D
            The simulation object
        """
        self.ax_energy.clear()

        # Get energy data
        times = np.arange(len(sim.energies['total'])) * sim.dt

        if len(times) > 0:
            # Plot different energy components
            self.ax_energy.plot(times, sim.energies['kinetic'],
                              'b-', label='Kinetic', linewidth=1.5)
            self.ax_energy.plot(times, sim.energies['potential'],
                              'r-', label='Potential', linewidth=1.5)
            self.ax_energy.plot(times, sim.energies['total'],
                              'k-', label='Total', linewidth=2)

            if max(sim.energies['bond']) > 0.01:
                self.ax_energy.plot(times, sim.energies['bond'],
                                  'g--', label='Bond', alpha=0.7)

            if max(sim.energies['lj']) != 0:
                self.ax_energy.plot(times, sim.energies['lj'],
                                  'm--', label='LJ', alpha=0.7)

            if max(sim.energies['angle']) > 0.01:
                self.ax_energy.plot(times, sim.energies['angle'],
                                  'c--', label='Angle', alpha=0.7)

        self.ax_energy.set_xlabel('Time (ps)')
        self.ax_energy.set_ylabel('Energy (kJ/mol)')
        self.ax_energy.set_title('Energy Evolution')
        self.ax_energy.legend(loc='best', fontsize=8)
        self.ax_energy.grid(True, alpha=0.3)

    def render(self, sim, show_bonds: bool = True):
        """
        Complete rendering: simulation view + energy plot.

        Parameters
        ----------
        sim : Simulation2D
            The simulation object to render
        show_bonds : bool
            Whether to display bonds
        """
        # Update trails
        if self.show_trails:
            positions = sim.get_positions()
            if len(self.trails) == 0:
                self.trails = [[] for _ in range(len(positions))]

            for i, pos in enumerate(positions):
                if i < len(self.trails):
                    self.trails[i].append(pos.copy())
                    if len(self.trails[i]) > self.trail_length:
                        self.trails[i].pop(0)

        # Render simulation frame
        self.render_frame(sim, show_bonds=show_bonds)

        # Update energy plot
        self.update_energy_plot(sim)

        plt.tight_layout()
        plt.draw()
        plt.pause(0.001)

    def save_frame(self, filename: str):
        """Save current frame to file."""
        self.fig.savefig(filename, dpi=150, bbox_inches='tight')

    def create_animation(self, sim, n_steps: int, interval: int = 50,
                        filename: Optional[str] = None):
        """
        Create an animation of the simulation.

        Parameters
        ----------
        sim : Simulation2D
            The simulation object
        n_steps : int
            Number of steps to animate
        interval : int
            Delay between frames in milliseconds
        filename : str, optional
            If provided, save animation to file
        """
        def animate(frame):
            # Run simulation step
            sim.integrate_step()

            # Render
            self.render(sim, show_bonds=True)

            return self.ax_sim.artists + self.ax_energy.lines

        anim = animation.FuncAnimation(self.fig, animate, frames=n_steps,
                                      interval=interval, blit=False)

        if filename is not None:
            if filename.endswith('.gif'):
                anim.save(filename, writer='pillow', fps=20)
            elif filename.endswith('.mp4'):
                anim.save(filename, writer='ffmpeg', fps=20)

        return anim

    def show(self):
        """Display the plot window."""
        plt.show()

    def close(self):
        """Close the plot window."""
        plt.close(self.fig)


def plot_energy_summary(sim, figsize: Tuple[float, float] = (12, 8)):
    """
    Create a comprehensive energy analysis plot.

    Parameters
    ----------
    sim : Simulation2D
        The simulation object
    figsize : Tuple[float, float]
        Figure size in inches
    """
    fig, axes = plt.subplots(2, 2, figsize=figsize)

    times = np.arange(len(sim.energies['total'])) * sim.dt

    # Total energy
    axes[0, 0].plot(times, sim.energies['total'], 'k-', linewidth=2)
    axes[0, 0].set_xlabel('Time (ps)')
    axes[0, 0].set_ylabel('Total Energy (kJ/mol)')
    axes[0, 0].set_title('Energy Conservation')
    axes[0, 0].grid(True, alpha=0.3)

    # Kinetic vs Potential
    axes[0, 1].plot(times, sim.energies['kinetic'], 'b-', label='Kinetic', linewidth=1.5)
    axes[0, 1].plot(times, sim.energies['potential'], 'r-', label='Potential', linewidth=1.5)
    axes[0, 1].set_xlabel('Time (ps)')
    axes[0, 1].set_ylabel('Energy (kJ/mol)')
    axes[0, 1].set_title('Kinetic vs Potential Energy')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)

    # Components
    axes[1, 0].plot(times, sim.energies['lj'], 'm-', label='Lennard-Jones', linewidth=1.5)
    axes[1, 0].plot(times, sim.energies['bond'], 'g-', label='Bonds', linewidth=1.5)
    axes[1, 0].plot(times, sim.energies['angle'], 'c-', label='Angles', linewidth=1.5)
    axes[1, 0].set_xlabel('Time (ps)')
    axes[1, 0].set_ylabel('Energy (kJ/mol)')
    axes[1, 0].set_title('Energy Components')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)

    # Temperature
    if sim.energies['kinetic']:
        k_B = 0.00831446261815324
        n_df = 2 * len(sim.particles)
        temperatures = [2 * ke / (k_B * n_df) if n_df > 0 else 0
                       for ke in sim.energies['kinetic']]

        axes[1, 1].plot(times, temperatures, 'orange', linewidth=2)
        if sim.target_temperature is not None:
            axes[1, 1].axhline(sim.target_temperature, color='red',
                             linestyle='--', label=f'Target: {sim.target_temperature}K')
            axes[1, 1].legend()
        axes[1, 1].set_xlabel('Time (ps)')
        axes[1, 1].set_ylabel('Temperature (K)')
        axes[1, 1].set_title('Temperature Evolution')
        axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    return fig
