"""
Educational 2D/3D Simulations using GROMOS-RS Interactions
===========================================================

This package provides educational molecular dynamics simulations
using the refined interactions from gromos-rs.

Modules:
- simulation2d: Core 2D simulation engine
- simulation3d: Core 3D simulation engine
- renderer: Real-time 2D matplotlib visualization
- renderer3d: Interactive 3D matplotlib visualization
- interactions: 2D/3D interaction wrappers for gromos-rs
- movie_utils: Movie and animation creation utilities
"""

# 2D modules
from .simulation2d import Simulation2D, Particle
from .renderer import Renderer2D
from .interactions import LennardJones2D, HarmonicBond2D, AnglePotential2D

# 3D modules
from .simulation3d import Simulation3D, Particle3D
from .renderer3d import Renderer3D, InteractiveRenderer3D

# Movie utilities
from .movie_utils import (
    MovieRecorder,
    record_simulation_movie,
    make_quick_movie,
    make_high_quality_movie,
    make_gif_animation,
    make_rotating_3d_movie
)

__all__ = [
    # 2D
    'Simulation2D',
    'Particle',
    'Renderer2D',
    'LennardJones2D',
    'HarmonicBond2D',
    'AnglePotential2D',
    # 3D
    'Simulation3D',
    'Particle3D',
    'Renderer3D',
    'InteractiveRenderer3D',
    # Movies
    'MovieRecorder',
    'record_simulation_movie',
    'make_quick_movie',
    'make_high_quality_movie',
    'make_gif_animation',
    'make_rotating_3d_movie',
]

__version__ = '2.0.0'
