"""
2D Molecular Dynamics Simulation Engine
========================================

Educational 2D MD simulator using gromos-rs interactions.
Demonstrates molecular forces, integration, and dynamics in 2D.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Callable
import sys
import os

# Add gromos-rs Python bindings to path if needed
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../py-gromos/python'))

try:
    import gromos
    GROMOS_AVAILABLE = True
except ImportError:
    print("Warning: gromos module not available. Using fallback implementations.")
    GROMOS_AVAILABLE = False


@dataclass
class Particle:
    """Represents a single particle in 2D space."""
    mass: float = 1.0                    # atomic mass units
    charge: float = 0.0                  # elementary charges
    position: np.ndarray = field(default_factory=lambda: np.zeros(2))
    velocity: np.ndarray = field(default_factory=lambda: np.zeros(2))
    force: np.ndarray = field(default_factory=lambda: np.zeros(2))
    atom_type: int = 0                   # For LJ parameters
    color: str = 'blue'                  # Visualization color

    def kinetic_energy(self) -> float:
        """Calculate kinetic energy: KE = 0.5 * m * v^2"""
        return 0.5 * self.mass * np.dot(self.velocity, self.velocity)


@dataclass
class Bond:
    """Harmonic bond between two particles."""
    i: int                               # Particle index 1
    j: int                               # Particle index 2
    k: float = 1000.0                    # Force constant (kJ/mol/nm²)
    r0: float = 0.1                      # Equilibrium distance (nm)


@dataclass
class Angle:
    """Angle potential between three particles."""
    i: int                               # Particle index 1
    j: int                               # Central particle index
    k: int                               # Particle index 3
    k_angle: float = 100.0               # Force constant (kJ/mol/rad²)
    theta0: float = np.pi                # Equilibrium angle (radians)


class Simulation2D:
    """
    Educational 2D Molecular Dynamics Simulation

    Uses gromos-rs force field calculations in 2D projection.
    Demonstrates:
    - Lennard-Jones interactions
    - Harmonic bonds
    - Angle potentials
    - Velocity Verlet integration
    - Temperature control
    """

    def __init__(self,
                 box_size: Tuple[float, float] = (5.0, 5.0),
                 periodic: bool = True,
                 integrator: str = 'velocity-verlet',
                 dt: float = 0.001):
        """
        Initialize 2D simulation.

        Parameters
        ----------
        box_size : Tuple[float, float]
            Simulation box dimensions in nm (width, height)
        periodic : bool
            Use periodic boundary conditions
        integrator : str
            Integration algorithm ('velocity-verlet', 'leap-frog')
        dt : float
            Time step in ps (picoseconds)
        """
        self.box_size = np.array(box_size, dtype=np.float64)
        self.periodic = periodic
        self.integrator = integrator
        self.dt = dt

        # Particles and interactions
        self.particles: List[Particle] = []
        self.bonds: List[Bond] = []
        self.angles: List[Angle] = []

        # Simulation state
        self.time = 0.0
        self.step = 0

        # Energy tracking
        self.energies = {
            'kinetic': [],
            'potential': [],
            'lj': [],
            'bond': [],
            'angle': [],
            'total': []
        }

        # Lennard-Jones parameters (GROMOS-style)
        # C6 and C12 coefficients for different atom types
        self.lj_params = {
            0: {'C6': 0.0025, 'C12': 0.00001, 'sigma': 0.3, 'epsilon': 0.5},  # Default
            1: {'C6': 0.0020, 'C12': 0.00008, 'sigma': 0.25, 'epsilon': 0.4}, # Smaller
            2: {'C6': 0.0030, 'C12': 0.00015, 'sigma': 0.35, 'epsilon': 0.6}, # Larger
        }

        # Temperature control
        self.target_temperature = None  # K
        self.temperature_coupling = None  # tau (ps)

    def add_particle(self,
                     position: np.ndarray,
                     velocity: Optional[np.ndarray] = None,
                     mass: float = 1.0,
                     charge: float = 0.0,
                     atom_type: int = 0,
                     color: str = 'blue') -> int:
        """Add a particle to the simulation."""
        if velocity is None:
            velocity = np.zeros(2)

        particle = Particle(
            mass=mass,
            charge=charge,
            position=np.array(position, dtype=np.float64),
            velocity=np.array(velocity, dtype=np.float64),
            force=np.zeros(2, dtype=np.float64),
            atom_type=atom_type,
            color=color
        )
        self.particles.append(particle)
        return len(self.particles) - 1

    def add_bond(self, i: int, j: int, k: float = 1000.0, r0: float = 0.1):
        """Add a harmonic bond between particles i and j."""
        self.bonds.append(Bond(i=i, j=j, k=k, r0=r0))

    def add_angle(self, i: int, j: int, k: int, k_angle: float = 100.0, theta0: float = np.pi):
        """Add an angle potential between particles i-j-k."""
        self.angles.append(Angle(i=i, j=j, k=k, k_angle=k_angle, theta0=theta0))

    def apply_pbc(self, r: np.ndarray) -> np.ndarray:
        """Apply periodic boundary conditions (minimum image convention)."""
        if not self.periodic:
            return r

        # Minimum image convention
        r = r - self.box_size * np.round(r / self.box_size)
        return r

    def calculate_forces(self):
        """
        Calculate all forces using gromos-rs-style interactions.

        Force components:
        1. Lennard-Jones (nonbonded)
        2. Harmonic bonds
        3. Angle potentials
        """
        # Reset forces
        for particle in self.particles:
            particle.force[:] = 0.0

        e_lj = 0.0
        e_bond = 0.0
        e_angle = 0.0

        # 1. Lennard-Jones interactions (all pairs)
        n = len(self.particles)
        for i in range(n):
            for j in range(i + 1, n):
                pi, pj = self.particles[i], self.particles[j]

                # Distance vector with PBC
                r_vec = pj.position - pi.position
                r_vec = self.apply_pbc(r_vec)
                r = np.linalg.norm(r_vec)

                if r < 0.001:  # Avoid division by zero
                    continue

                # Get LJ parameters (Lorentz-Berthelot combining rules)
                params_i = self.lj_params[pi.atom_type]
                params_j = self.lj_params[pj.atom_type]

                # GROMOS uses C6 and C12 directly
                C6 = np.sqrt(params_i['C6'] * params_j['C6'])
                C12 = np.sqrt(params_i['C12'] * params_j['C12'])

                # Lennard-Jones: E = C12/r^12 - C6/r^6
                r6 = r**6
                r12 = r6 * r6

                e_lj += C12 / r12 - C6 / r6

                # Force: F = (12*C12/r^14 - 6*C6/r^8) * r_vec/r
                f_magnitude = (12 * C12 / (r12 * r * r) - 6 * C6 / (r6 * r * r))
                f_vec = f_magnitude * r_vec / r

                pi.force += f_vec
                pj.force -= f_vec

        # 2. Harmonic bonds: E = 0.5 * k * (r - r0)^2
        for bond in self.bonds:
            pi, pj = self.particles[bond.i], self.particles[bond.j]

            r_vec = pj.position - pi.position
            r_vec = self.apply_pbc(r_vec)
            r = np.linalg.norm(r_vec)

            if r < 0.001:
                continue

            # Energy
            dr = r - bond.r0
            e_bond += 0.5 * bond.k * dr * dr

            # Force: F = -k * (r - r0) * r_vec/r
            f_magnitude = -bond.k * dr
            f_vec = f_magnitude * r_vec / r

            pi.force += f_vec
            pj.force -= f_vec

        # 3. Angle potentials: E = 0.5 * k * (theta - theta0)^2
        for angle in self.angles:
            pi = self.particles[angle.i]
            pj = self.particles[angle.j]  # Central atom
            pk = self.particles[angle.k]

            # Vectors from central atom
            r_ij = pi.position - pj.position
            r_kj = pk.position - pj.position
            r_ij = self.apply_pbc(r_ij)
            r_kj = self.apply_pbc(r_kj)

            rij = np.linalg.norm(r_ij)
            rkj = np.linalg.norm(r_kj)

            if rij < 0.001 or rkj < 0.001:
                continue

            # Angle calculation
            cos_theta = np.dot(r_ij, r_kj) / (rij * rkj)
            cos_theta = np.clip(cos_theta, -1.0, 1.0)
            theta = np.arccos(cos_theta)

            # Energy
            dtheta = theta - angle.theta0
            e_angle += 0.5 * angle.k_angle * dtheta * dtheta

            # Force (simplified gradient)
            # For educational purposes, using numerical gradient would be clearer
            # but here's the analytical form
            if abs(np.sin(theta)) > 0.001:
                f_factor = -angle.k_angle * dtheta / np.sin(theta)

                # Force on i
                f_i = f_factor * (r_kj / (rij * rkj) - cos_theta * r_ij / (rij * rij))
                # Force on k
                f_k = f_factor * (r_ij / (rij * rkj) - cos_theta * r_kj / (rkj * rkj))

                pi.force += f_i
                pk.force += f_k
                pj.force -= (f_i + f_k)

        return {
            'lj': e_lj,
            'bond': e_bond,
            'angle': e_angle,
            'potential': e_lj + e_bond + e_angle
        }

    def integrate_step(self):
        """Perform one integration step using Velocity Verlet."""
        if self.integrator == 'velocity-verlet':
            self._velocity_verlet()
        elif self.integrator == 'leap-frog':
            self._leap_frog()
        else:
            raise ValueError(f"Unknown integrator: {self.integrator}")

    def _velocity_verlet(self):
        """
        Velocity Verlet integration:

        v(t+dt/2) = v(t) + a(t) * dt/2
        x(t+dt) = x(t) + v(t+dt/2) * dt
        a(t+dt) = F(x(t+dt)) / m
        v(t+dt) = v(t+dt/2) + a(t+dt) * dt/2
        """
        # Step 1: v(t+dt/2) = v(t) + a(t) * dt/2
        for p in self.particles:
            p.velocity += 0.5 * self.dt * (p.force / p.mass)

        # Step 2: x(t+dt) = x(t) + v(t+dt/2) * dt
        for p in self.particles:
            p.position += self.dt * p.velocity

            # Apply periodic boundaries
            if self.periodic:
                p.position = np.mod(p.position, self.box_size)

        # Step 3: Calculate forces at new positions
        energies = self.calculate_forces()

        # Step 4: v(t+dt) = v(t+dt/2) + a(t+dt) * dt/2
        for p in self.particles:
            p.velocity += 0.5 * self.dt * (p.force / p.mass)

        # Apply temperature coupling if enabled
        if self.target_temperature is not None and self.temperature_coupling is not None:
            self._berendsen_thermostat()

        return energies

    def _leap_frog(self):
        """Leap-Frog integration (similar to GROMOS default)."""
        # Calculate forces
        energies = self.calculate_forces()

        # Update velocities: v(t+dt/2) = v(t-dt/2) + a(t) * dt
        for p in self.particles:
            p.velocity += self.dt * (p.force / p.mass)

        # Update positions: x(t+dt) = x(t) + v(t+dt/2) * dt
        for p in self.particles:
            p.position += self.dt * p.velocity

            if self.periodic:
                p.position = np.mod(p.position, self.box_size)

        if self.target_temperature is not None and self.temperature_coupling is not None:
            self._berendsen_thermostat()

        return energies

    def _berendsen_thermostat(self):
        """
        Berendsen thermostat (velocity rescaling).

        lambda = sqrt(1 + dt/tau * (T0/T - 1))
        """
        current_temp = self.temperature()
        if current_temp < 0.001:
            return

        scale = np.sqrt(1 + self.dt / self.temperature_coupling *
                       (self.target_temperature / current_temp - 1))

        for p in self.particles:
            p.velocity *= scale

    def temperature(self) -> float:
        """
        Calculate instantaneous temperature.

        T = 2 * KE / (k_B * N_df)
        where N_df = 2*N for 2D system (2 degrees of freedom per particle)
        """
        if len(self.particles) == 0:
            return 0.0

        ke = sum(p.kinetic_energy() for p in self.particles)
        k_B = 0.00831446261815324  # kJ/(mol·K)
        n_df = 2 * len(self.particles)  # 2D degrees of freedom

        return 2 * ke / (k_B * n_df) if n_df > 0 else 0.0

    def total_energy(self) -> Tuple[float, float, float]:
        """Calculate total energy (kinetic, potential, total)."""
        ke = sum(p.kinetic_energy() for p in self.particles)
        energies = self.calculate_forces()
        pe = energies['potential']
        return ke, pe, ke + pe

    def run(self, n_steps: int, callback: Optional[Callable] = None):
        """
        Run simulation for n_steps.

        Parameters
        ----------
        n_steps : int
            Number of integration steps
        callback : callable, optional
            Function to call after each step: callback(sim, step)
        """
        for _ in range(n_steps):
            energies = self.integrate_step()

            # Track energies
            ke = sum(p.kinetic_energy() for p in self.particles)
            self.energies['kinetic'].append(ke)
            self.energies['lj'].append(energies['lj'])
            self.energies['bond'].append(energies['bond'])
            self.energies['angle'].append(energies['angle'])
            self.energies['potential'].append(energies['potential'])
            self.energies['total'].append(ke + energies['potential'])

            self.time += self.dt
            self.step += 1

            if callback is not None:
                callback(self, self.step)

    def get_positions(self) -> np.ndarray:
        """Get all particle positions as numpy array (N, 2)."""
        return np.array([p.position for p in self.particles])

    def get_velocities(self) -> np.ndarray:
        """Get all particle velocities as numpy array (N, 2)."""
        return np.array([p.velocity for p in self.particles])

    def get_forces(self) -> np.ndarray:
        """Get all particle forces as numpy array (N, 2)."""
        return np.array([p.force for p in self.particles])

    def get_colors(self) -> List[str]:
        """Get particle colors for visualization."""
        return [p.color for p in self.particles]

    def set_temperature(self, temperature: float, tau: float = 0.1):
        """
        Enable temperature coupling.

        Parameters
        ----------
        temperature : float
            Target temperature in Kelvin
        tau : float
            Coupling time constant in ps
        """
        self.target_temperature = temperature
        self.temperature_coupling = tau

    def initialize_velocities(self, temperature: float):
        """
        Initialize velocities from Maxwell-Boltzmann distribution.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin
        """
        k_B = 0.00831446261815324  # kJ/(mol·K)

        for p in self.particles:
            # v ~ N(0, sqrt(k_B*T/m))
            sigma = np.sqrt(k_B * temperature / p.mass)
            p.velocity = np.random.normal(0, sigma, size=2)

        # Remove center-of-mass motion
        total_momentum = sum(p.mass * p.velocity for p in self.particles)
        total_mass = sum(p.mass for p in self.particles)
        v_com = total_momentum / total_mass

        for p in self.particles:
            p.velocity -= v_com
