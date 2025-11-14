"""
3D Visualization and Rendering
===============================

Interactive 3D matplotlib visualization for molecular dynamics simulations.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from typing import Optional, Tuple


class Renderer3D:
    """
    Real-time 3D visualization for MD simulations.

    Features:
    - Interactive 3D particle rendering
    - Bond visualization
    - Rotatable view
    - Energy plots (2D subplot)
    - Multiple color schemes
    """

    def __init__(self,
                 box_size: Tuple[float, float, float],
                 figsize: Tuple[float, float] = (14, 6),
                 azimuth: float = 45,
                 elevation: float = 30):
        """
        Initialize 3D renderer.

        Parameters
        ----------
        box_size : Tuple[float, float, float]
            Simulation box size (x, y, z) in nm
        figsize : Tuple[float, float]
            Figure size in inches
        azimuth : float
            Initial azimuth viewing angle (degrees)
        elevation : float
            Initial elevation viewing angle (degrees)
        """
        self.box_size = np.array(box_size)
        self.figsize = figsize
        self.azimuth = azimuth
        self.elevation = elevation

        # Create figure with 3D subplot and energy subplot
        self.fig = plt.figure(figsize=figsize)
        self.ax_3d = self.fig.add_subplot(121, projection='3d')
        self.ax_energy = self.fig.add_subplot(122)

        # Setup 3D view
        self._setup_3d_axes()

        # Setup energy plot
        self.ax_energy.set_xlabel('Time (ps)')
        self.ax_energy.set_ylabel('Energy (kJ/mol)')
        self.ax_energy.set_title('Energy Evolution')
        self.ax_energy.grid(True, alpha=0.3)

        plt.tight_layout()

    def _setup_3d_axes(self):
        """Setup 3D axes with box."""
        ax = self.ax_3d

        # Set limits
        ax.set_xlim(0, self.box_size[0])
        ax.set_ylim(0, self.box_size[1])
        ax.set_zlim(0, self.box_size[2])

        # Labels
        ax.set_xlabel('x (nm)')
        ax.set_ylabel('y (nm)')
        ax.set_zlabel('z (nm)')

        # Set view angle
        ax.view_init(elev=self.elevation, azim=self.azimuth)

        # Equal aspect ratio
        ax.set_box_aspect([1, 1, 1])

        # Draw box edges
        self._draw_box()

    def _draw_box(self):
        """Draw simulation box edges."""
        ax = self.ax_3d
        x, y, z = self.box_size

        # Define box edges
        edges = [
            # Bottom face
            [(0, 0, 0), (x, 0, 0)],
            [(x, 0, 0), (x, y, 0)],
            [(x, y, 0), (0, y, 0)],
            [(0, y, 0), (0, 0, 0)],
            # Top face
            [(0, 0, z), (x, 0, z)],
            [(x, 0, z), (x, y, z)],
            [(x, y, z), (0, y, z)],
            [(0, y, z), (0, 0, z)],
            # Vertical edges
            [(0, 0, 0), (0, 0, z)],
            [(x, 0, 0), (x, 0, z)],
            [(x, y, 0), (x, y, z)],
            [(0, y, 0), (0, y, z)],
        ]

        for edge in edges:
            points = np.array(edge)
            ax.plot3D(points[:, 0], points[:, 1], points[:, 2],
                     'k-', linewidth=1, alpha=0.3)

    def render_frame(self, sim, show_bonds: bool = True):
        """
        Render a single frame of the 3D simulation.

        Parameters
        ----------
        sim : Simulation3D
            The simulation object to render
        show_bonds : bool
            Whether to display bonds
        """
        positions = sim.get_positions()
        colors = sim.get_colors()

        # Clear previous frame
        self.ax_3d.clear()
        self._setup_3d_axes()

        # Draw bonds
        if show_bonds and len(sim.bonds) > 0:
            for bond in sim.bonds:
                pi = positions[bond.i]
                pj = positions[bond.j]

                # Handle periodic boundaries
                dr = pj - pi
                if sim.periodic:
                    dr = dr - self.box_size * np.round(dr / self.box_size)
                    pj_wrapped = pi + dr

                    self.ax_3d.plot3D([pi[0], pj_wrapped[0]],
                                     [pi[1], pj_wrapped[1]],
                                     [pi[2], pj_wrapped[2]],
                                     'gray', linewidth=2, alpha=0.6)
                else:
                    self.ax_3d.plot3D([pi[0], pj[0]],
                                     [pi[1], pj[1]],
                                     [pi[2], pj[2]],
                                     'gray', linewidth=2, alpha=0.6)

        # Draw particles
        self.ax_3d.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
                          c=colors, s=100, alpha=0.8, edgecolors='black',
                          linewidth=1)

        # Display info
        info_text = f"Step: {sim.step}\n"
        info_text += f"Time: {sim.time:.3f} ps\n"
        info_text += f"T: {sim.temperature():.1f} K\n"
        info_text += f"Particles: {len(sim.particles)}"

        self.ax_3d.text2D(0.02, 0.98, info_text,
                         transform=self.ax_3d.transAxes,
                         verticalalignment='top',
                         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                         fontsize=9)

        self.ax_3d.set_title('3D Molecular Dynamics')

    def update_energy_plot(self, sim):
        """Update the energy evolution plot."""
        self.ax_energy.clear()

        times = np.arange(len(sim.energies['total'])) * sim.dt

        if len(times) > 0:
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

        self.ax_energy.set_xlabel('Time (ps)')
        self.ax_energy.set_ylabel('Energy (kJ/mol)')
        self.ax_energy.set_title('Energy Evolution')
        self.ax_energy.legend(loc='best', fontsize=8)
        self.ax_energy.grid(True, alpha=0.3)

    def render(self, sim, show_bonds: bool = True, azimuth: Optional[float] = None,
              elevation: Optional[float] = None):
        """
        Complete rendering: 3D view + energy plot.

        Parameters
        ----------
        sim : Simulation3D
            The simulation object to render
        show_bonds : bool
            Whether to display bonds
        azimuth : float, optional
            Viewing azimuth angle (if provided, updates view)
        elevation : float, optional
            Viewing elevation angle (if provided, updates view)
        """
        # Update view angles if provided
        if azimuth is not None:
            self.azimuth = azimuth
        if elevation is not None:
            self.elevation = elevation

        self.ax_3d.view_init(elev=self.elevation, azim=self.azimuth)

        # Render
        self.render_frame(sim, show_bonds=show_bonds)
        self.update_energy_plot(sim)

        plt.tight_layout()
        plt.draw()
        plt.pause(0.001)

    def save_frame(self, filename: str):
        """Save current frame to file."""
        self.fig.savefig(filename, dpi=150, bbox_inches='tight')

    def show(self):
        """Display the plot window."""
        plt.show()

    def close(self):
        """Close the plot window."""
        plt.close(self.fig)


class InteractiveRenderer3D:
    """
    Interactive 3D renderer with mouse controls.

    Features:
    - Click and drag to rotate
    - Scroll to zoom
    - Real-time updates
    """

    def __init__(self, box_size: Tuple[float, float, float], figsize: Tuple[float, float] = (10, 8)):
        """
        Initialize interactive 3D renderer.

        Parameters
        ----------
        box_size : Tuple[float, float, float]
            Simulation box size
        figsize : Tuple[float, float]
            Figure size
        """
        self.box_size = np.array(box_size)

        # Create figure
        self.fig = plt.figure(figsize=figsize)
        self.ax = self.fig.add_subplot(111, projection='3d')

        # Setup axes
        self.ax.set_xlim(0, box_size[0])
        self.ax.set_ylim(0, box_size[1])
        self.ax.set_zlim(0, box_size[2])
        self.ax.set_xlabel('x (nm)')
        self.ax.set_ylabel('y (nm)')
        self.ax.set_zlabel('z (nm)')
        self.ax.set_box_aspect([1, 1, 1])

        # Draw box
        self._draw_box()

        # Enable mouse interaction
        self.fig.canvas.mpl_connect('button_press_event', self._on_click)
        self.fig.canvas.mpl_connect('motion_notify_event', self._on_motion)
        self.fig.canvas.mpl_connect('scroll_event', self._on_scroll)

        self.dragging = False
        self.last_mouse = None

    def _draw_box(self):
        """Draw simulation box."""
        x, y, z = self.box_size

        edges = [
            [(0, 0, 0), (x, 0, 0)], [(x, 0, 0), (x, y, 0)],
            [(x, y, 0), (0, y, 0)], [(0, y, 0), (0, 0, 0)],
            [(0, 0, z), (x, 0, z)], [(x, 0, z), (x, y, z)],
            [(x, y, z), (0, y, z)], [(0, y, z), (0, 0, z)],
            [(0, 0, 0), (0, 0, z)], [(x, 0, 0), (x, 0, z)],
            [(x, y, 0), (x, y, z)], [(0, y, 0), (0, y, z)],
        ]

        for edge in edges:
            points = np.array(edge)
            self.ax.plot3D(points[:, 0], points[:, 1], points[:, 2],
                          'k-', linewidth=1, alpha=0.3)

    def _on_click(self, event):
        """Handle mouse click."""
        if event.inaxes == self.ax:
            self.dragging = True
            self.last_mouse = (event.xdata, event.ydata)

    def _on_motion(self, event):
        """Handle mouse motion (rotation)."""
        if self.dragging and event.inaxes == self.ax:
            if self.last_mouse is not None:
                # Rotate view
                dx = event.xdata - self.last_mouse[0]
                dy = event.ydata - self.last_mouse[1]

                # Update view angles
                azim, elev = self.ax.azim, self.ax.elev
                self.ax.view_init(elev=elev + dy * 0.5, azim=azim - dx * 0.5)

                self.fig.canvas.draw_idle()
                self.last_mouse = (event.xdata, event.ydata)

    def _on_scroll(self, event):
        """Handle scroll (zoom)."""
        if event.inaxes == self.ax:
            # Zoom in/out
            scale_factor = 1.1 if event.button == 'up' else 0.9

            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            zlim = self.ax.get_zlim()

            center_x = (xlim[0] + xlim[1]) / 2
            center_y = (ylim[0] + ylim[1]) / 2
            center_z = (zlim[0] + zlim[1]) / 2

            width_x = (xlim[1] - xlim[0]) * scale_factor
            width_y = (ylim[1] - ylim[0]) * scale_factor
            width_z = (zlim[1] - zlim[0]) * scale_factor

            self.ax.set_xlim(center_x - width_x/2, center_x + width_x/2)
            self.ax.set_ylim(center_y - width_y/2, center_y + width_y/2)
            self.ax.set_zlim(center_z - width_z/2, center_z + width_z/2)

            self.fig.canvas.draw_idle()

    def render(self, sim, show_bonds: bool = True):
        """Render simulation frame."""
        positions = sim.get_positions()
        colors = sim.get_colors()

        self.ax.clear()
        self._draw_box()

        # Restore limits
        self.ax.set_xlim(0, self.box_size[0])
        self.ax.set_ylim(0, self.box_size[1])
        self.ax.set_zlim(0, self.box_size[2])
        self.ax.set_xlabel('x (nm)')
        self.ax.set_ylabel('y (nm)')
        self.ax.set_zlabel('z (nm)')

        # Draw bonds
        if show_bonds:
            for bond in sim.bonds:
                pi = positions[bond.i]
                pj = positions[bond.j]
                self.ax.plot3D([pi[0], pj[0]], [pi[1], pj[1]], [pi[2], pj[2]],
                              'gray', linewidth=2, alpha=0.6)

        # Draw particles
        self.ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
                       c=colors, s=100, alpha=0.8, edgecolors='black')

        plt.draw()
        plt.pause(0.001)

    def show(self):
        """Show the interactive plot."""
        plt.show()
