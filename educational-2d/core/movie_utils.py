"""
Movie Rendering Utilities
==========================

Easy movie creation from 2D/3D molecular dynamics simulations.
Supports MP4, GIF, and frame sequence export.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from typing import Optional, Callable, List
import os
from pathlib import Path


class MovieRecorder:
    """
    Records simulation frames and exports to video formats.

    Supports:
    - MP4 video (requires ffmpeg)
    - Animated GIF (requires pillow)
    - Frame sequence (PNG/JPG)
    """

    def __init__(self,
                 output_file: str = "simulation.mp4",
                 fps: int = 30,
                 dpi: int = 100,
                 bitrate: int = 1800):
        """
        Initialize movie recorder.

        Parameters
        ----------
        output_file : str
            Output filename (.mp4, .gif, or directory for frames)
        fps : int
            Frames per second
        dpi : int
            Resolution (dots per inch)
        bitrate : int
            Video bitrate for MP4 (higher = better quality)
        """
        self.output_file = output_file
        self.fps = fps
        self.dpi = dpi
        self.bitrate = bitrate

        # Determine output format
        self.format = self._detect_format()

        # Frame storage
        self.frames = []
        self.frame_count = 0

    def _detect_format(self) -> str:
        """Detect output format from filename."""
        if self.output_file.endswith('.mp4'):
            return 'mp4'
        elif self.output_file.endswith('.gif'):
            return 'gif'
        elif os.path.isdir(self.output_file) or '/' in self.output_file:
            return 'frames'
        else:
            # Default to MP4
            self.output_file += '.mp4'
            return 'mp4'

    def add_frame(self, fig):
        """
        Add a matplotlib figure as a frame.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            Figure to capture
        """
        # Convert figure to image array
        fig.canvas.draw()
        frame = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
        frame = frame.reshape(fig.canvas.get_width_height()[::-1] + (3,))

        self.frames.append(frame)
        self.frame_count += 1

    def save(self):
        """Save recorded frames to file."""
        if len(self.frames) == 0:
            print("Warning: No frames to save!")
            return

        print(f"\nSaving {self.frame_count} frames to {self.output_file}...")

        if self.format == 'mp4':
            self._save_mp4()
        elif self.format == 'gif':
            self._save_gif()
        elif self.format == 'frames':
            self._save_frames()

    def _save_mp4(self):
        """Save as MP4 video using matplotlib animation."""
        try:
            # Create figure for animation
            fig = plt.figure(figsize=(self.frames[0].shape[1]/self.dpi,
                                     self.frames[0].shape[0]/self.dpi))
            ax = fig.add_subplot(111)
            ax.axis('off')

            # Display first frame
            im = ax.imshow(self.frames[0])

            def update(frame_idx):
                im.set_data(self.frames[frame_idx])
                return [im]

            # Create animation
            anim = animation.FuncAnimation(
                fig, update, frames=len(self.frames),
                interval=1000/self.fps, blit=True
            )

            # Save
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=self.fps, bitrate=self.bitrate)
            anim.save(self.output_file, writer=writer, dpi=self.dpi)

            plt.close(fig)
            print(f"✓ Saved MP4: {self.output_file}")
            print(f"  Frames: {self.frame_count}, FPS: {self.fps}, Duration: {self.frame_count/self.fps:.2f}s")

        except Exception as e:
            print(f"Error saving MP4: {e}")
            print("Tip: Install ffmpeg: conda install -c conda-forge ffmpeg")

    def _save_gif(self):
        """Save as animated GIF."""
        try:
            from PIL import Image

            # Convert frames to PIL Images
            images = [Image.fromarray(frame) for frame in self.frames]

            # Save as GIF
            images[0].save(
                self.output_file,
                save_all=True,
                append_images=images[1:],
                duration=int(1000/self.fps),
                loop=0,
                optimize=True
            )

            print(f"✓ Saved GIF: {self.output_file}")
            print(f"  Frames: {self.frame_count}, FPS: {self.fps}")

        except ImportError:
            print("Error: Pillow (PIL) not installed")
            print("Install with: pip install Pillow")
        except Exception as e:
            print(f"Error saving GIF: {e}")

    def _save_frames(self):
        """Save as sequence of PNG frames."""
        try:
            # Create output directory
            Path(self.output_file).mkdir(parents=True, exist_ok=True)

            # Save each frame
            from PIL import Image
            for i, frame in enumerate(self.frames):
                img = Image.fromarray(frame)
                filename = os.path.join(self.output_file, f"frame_{i:04d}.png")
                img.save(filename)

            print(f"✓ Saved {self.frame_count} frames to {self.output_file}/")

        except Exception as e:
            print(f"Error saving frames: {e}")


def record_simulation_movie(
    sim,
    renderer,
    n_steps: int,
    output_file: str = "simulation.mp4",
    render_interval: int = 1,
    fps: int = 30,
    dpi: int = 100,
    show_bonds: bool = True,
    progress_interval: int = 100
):
    """
    Record a simulation as a movie.

    Parameters
    ----------
    sim : Simulation2D or Simulation3D
        Simulation object
    renderer : Renderer2D or Renderer3D
        Renderer object
    n_steps : int
        Number of simulation steps to run
    output_file : str
        Output filename (.mp4 or .gif)
    render_interval : int
        Render every N steps (1 = every step)
    fps : int
        Frames per second in output video
    dpi : int
        Resolution
    show_bonds : bool
        Whether to show bonds (if applicable)
    progress_interval : int
        Print progress every N steps

    Returns
    -------
    MovieRecorder
        Recorder object (already saved)
    """
    recorder = MovieRecorder(output_file=output_file, fps=fps, dpi=dpi)

    print(f"\nRecording simulation movie...")
    print(f"  Steps: {n_steps}")
    print(f"  Render interval: {render_interval}")
    print(f"  Expected frames: {n_steps // render_interval}")
    print(f"  Output: {output_file}")
    print()

    for step in range(n_steps):
        # Run simulation step
        sim.integrate_step()

        # Render and record frame
        if step % render_interval == 0:
            renderer.render(sim, show_bonds=show_bonds)
            recorder.add_frame(renderer.fig)

        # Progress
        if step % progress_interval == 0:
            print(f"  Step {step}/{n_steps} ({100*step/n_steps:.1f}%) - "
                  f"Frames: {recorder.frame_count}")

    print(f"  Step {n_steps}/{n_steps} (100.0%) - Frames: {recorder.frame_count}")

    # Save movie
    recorder.save()

    return recorder


def create_comparison_movie(
    sims: List,
    renderers: List,
    titles: List[str],
    n_steps: int,
    output_file: str = "comparison.mp4",
    render_interval: int = 1,
    fps: int = 30,
    dpi: int = 100
):
    """
    Create a side-by-side comparison movie of multiple simulations.

    Parameters
    ----------
    sims : List
        List of simulation objects
    renderers : List
        List of renderer objects
    titles : List[str]
        Title for each simulation
    n_steps : int
        Number of steps
    output_file : str
        Output filename
    render_interval : int
        Render every N steps
    fps : int
        Frames per second
    dpi : int
        Resolution

    Returns
    -------
    MovieRecorder
        Recorder object
    """
    n_sims = len(sims)
    recorder = MovieRecorder(output_file=output_file, fps=fps, dpi=dpi)

    print(f"\nRecording comparison movie ({n_sims} simulations)...")

    # Create combined figure
    fig, axes = plt.subplots(1, n_sims, figsize=(6*n_sims, 6))
    if n_sims == 1:
        axes = [axes]

    for ax, title in zip(axes, titles):
        ax.set_title(title)

    for step in range(n_steps):
        # Run all simulations
        for sim in sims:
            sim.integrate_step()

        # Render
        if step % render_interval == 0:
            for i, (sim, renderer) in enumerate(zip(sims, renderers)):
                plt.sca(axes[i])
                renderer.render_frame(sim, show_bonds=True)

            recorder.add_frame(fig)

        if step % 100 == 0:
            print(f"  Step {step}/{n_steps} - Frames: {recorder.frame_count}")

    plt.close(fig)
    recorder.save()

    return recorder


def make_rotating_3d_movie(
    sim,
    renderer_3d,
    n_frames: int = 360,
    output_file: str = "rotation.mp4",
    elevation: float = 30,
    fps: int = 30
):
    """
    Create a rotating view movie of a 3D simulation snapshot.

    Parameters
    ----------
    sim : Simulation3D
        3D simulation object
    renderer_3d : Renderer3D
        3D renderer
    n_frames : int
        Number of rotation frames (360 = full rotation)
    output_file : str
        Output file
    elevation : float
        Viewing elevation angle
    fps : int
        Frames per second

    Returns
    -------
    MovieRecorder
        Recorder object
    """
    recorder = MovieRecorder(output_file=output_file, fps=fps)

    print(f"\nCreating rotating 3D view ({n_frames} frames)...")

    for frame in range(n_frames):
        azimuth = 360 * frame / n_frames

        renderer_3d.render(sim, azimuth=azimuth, elevation=elevation)
        recorder.add_frame(renderer_3d.fig)

        if frame % 30 == 0:
            print(f"  Frame {frame}/{n_frames} (azimuth: {azimuth:.1f}°)")

    recorder.save()
    return recorder


# Convenience functions for common movie types

def make_quick_movie(sim, renderer, n_steps: int, filename: str = "movie.mp4"):
    """Quick movie with default settings."""
    return record_simulation_movie(
        sim, renderer, n_steps,
        output_file=filename,
        render_interval=max(1, n_steps // 300),  # Cap at ~300 frames
        fps=30
    )


def make_high_quality_movie(sim, renderer, n_steps: int, filename: str = "movie_hq.mp4"):
    """High quality movie (every frame, high DPI)."""
    return record_simulation_movie(
        sim, renderer, n_steps,
        output_file=filename,
        render_interval=1,
        fps=60,
        dpi=150
    )


def make_gif_animation(sim, renderer, n_steps: int, filename: str = "animation.gif"):
    """Animated GIF (good for web/presentations)."""
    return record_simulation_movie(
        sim, renderer, n_steps,
        output_file=filename,
        render_interval=max(1, n_steps // 100),  # GIFs should be smaller
        fps=15
    )


# Example usage
if __name__ == "__main__":
    print("Movie Utilities for Educational MD Simulations")
    print("=" * 50)
    print("\nAvailable functions:")
    print("  - record_simulation_movie(): Record any simulation")
    print("  - make_quick_movie(): Fast movie with good quality")
    print("  - make_high_quality_movie(): High FPS/DPI movie")
    print("  - make_gif_animation(): Animated GIF")
    print("  - create_comparison_movie(): Side-by-side comparison")
    print("  - make_rotating_3d_movie(): Rotating 3D view")
    print("\nExample:")
    print("  from core.movie_utils import make_quick_movie")
    print("  make_quick_movie(sim, renderer, n_steps=1000, filename='my_movie.mp4')")
