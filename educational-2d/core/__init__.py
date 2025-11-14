"""
Educational 2D Simulations using GROMOS-RS Interactions
========================================================

This package provides educational 2D molecular dynamics simulations
using the refined interactions from gromos-rs.

Modules:
- simulation2d: Core 2D simulation engine
- renderer: Real-time matplotlib visualization
- interactions: 2D interaction wrappers for gromos-rs
"""

from .simulation2d import Simulation2D, Particle
from .renderer import Renderer2D
from .interactions import LennardJones2D, HarmonicBond2D, AnglePotential2D

__all__ = [
    'Simulation2D',
    'Particle',
    'Renderer2D',
    'LennardJones2D',
    'HarmonicBond2D',
    'AnglePotential2D',
]

__version__ = '1.0.0'
