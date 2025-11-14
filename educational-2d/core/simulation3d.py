"""
3D Molecular Dynamics Simulation Engine
========================================

Simple 3D MD simulator using gromos-rs style interactions.
Extension of the 2D framework to 3 dimensions.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Callable


@dataclass
class Particle3D:
    """Represents a single particle in 3D space."""
    mass: float = 1.0
    charge: float = 0.0
    position: np.ndarray = field(default_factory=lambda: np.zeros(3))
    velocity: np.ndarray = field(default_factory=lambda: np.zeros(3))
    force: np.ndarray = field(default_factory=lambda: np.zeros(3))
    atom_type: int = 0
    color: str = 'blue'

    def kinetic_energy(self) -> float:
        """Calculate kinetic energy: KE = 0.5 * m * v^2"""
        return 0.5 * self.mass * np.dot(self.velocity, self.velocity)


@dataclass
class Bond3D:
    """Harmonic bond between two particles."""
    i: int
    j: int
    k: float = 1000.0
    r0: float = 0.1


@dataclass
class Angle3D:
    """Angle potential between three particles."""
    i: int
    j: int
    k: int
    k_angle: float = 100.0
    theta0: float = np.pi


class Simulation3D:
    """
    Educational 3D Molecular Dynamics Simulation

    Similar to Simulation2D but in 3D space.
    """

    def __init__(self,
                 box_size: Tuple[float, float, float] = (5.0, 5.0, 5.0),
                 periodic: bool = True,
                 integrator: str = 'velocity-verlet',
                 dt: float = 0.001):
        """
        Initialize 3D simulation.

        Parameters
        ----------
        box_size : Tuple[float, float, float]
            Simulation box dimensions in nm (x, y, z)
        periodic : bool
            Use periodic boundary conditions
        integrator : str
            Integration algorithm
        dt : float
            Time step in ps
        """
        self.box_size = np.array(box_size, dtype=np.float64)
        self.periodic = periodic
        self.integrator = integrator
        self.dt = dt

        self.particles: List[Particle3D] = []
        self.bonds: List[Bond3D] = []
        self.angles: List[Angle3D] = []

        self.time = 0.0
        self.step = 0

        self.energies = {
            'kinetic': [],
            'potential': [],
            'lj': [],
            'bond': [],
            'angle': [],
            'total': []
        }

        # LJ parameters
        self.lj_params = {
            0: {'C6': 0.0025, 'C12': 0.00001, 'sigma': 0.3, 'epsilon': 0.5},
            1: {'C6': 0.0020, 'C12': 0.00008, 'sigma': 0.25, 'epsilon': 0.4},
            2: {'C6': 0.0030, 'C12': 0.00015, 'sigma': 0.35, 'epsilon': 0.6},
        }

        self.target_temperature = None
        self.temperature_coupling = None

    def add_particle(self,
                     position: np.ndarray,
                     velocity: Optional[np.ndarray] = None,
                     mass: float = 1.0,
                     charge: float = 0.0,
                     atom_type: int = 0,
                     color: str = 'blue') -> int:
        """Add a particle to the simulation."""
        if velocity is None:
            velocity = np.zeros(3)

        particle = Particle3D(
            mass=mass,
            charge=charge,
            position=np.array(position, dtype=np.float64),
            velocity=np.array(velocity, dtype=np.float64),
            force=np.zeros(3, dtype=np.float64),
            atom_type=atom_type,
            color=color
        )
        self.particles.append(particle)
        return len(self.particles) - 1

    def add_bond(self, i: int, j: int, k: float = 1000.0, r0: float = 0.1):
        """Add a harmonic bond."""
        self.bonds.append(Bond3D(i=i, j=j, k=k, r0=r0))

    def add_angle(self, i: int, j: int, k: int, k_angle: float = 100.0, theta0: float = np.pi):
        """Add an angle potential."""
        self.angles.append(Angle3D(i=i, j=j, k=k, k_angle=k_angle, theta0=theta0))

    def apply_pbc(self, r: np.ndarray) -> np.ndarray:
        """Apply periodic boundary conditions."""
        if not self.periodic:
            return r
        return r - self.box_size * np.round(r / self.box_size)

    def calculate_forces(self):
        """Calculate all forces."""
        for particle in self.particles:
            particle.force[:] = 0.0

        e_lj = 0.0
        e_bond = 0.0
        e_angle = 0.0

        # Lennard-Jones
        n = len(self.particles)
        for i in range(n):
            for j in range(i + 1, n):
                pi, pj = self.particles[i], self.particles[j]

                r_vec = pj.position - pi.position
                r_vec = self.apply_pbc(r_vec)
                r = np.linalg.norm(r_vec)

                if r < 0.001:
                    continue

                params_i = self.lj_params[pi.atom_type]
                params_j = self.lj_params[pj.atom_type]

                C6 = np.sqrt(params_i['C6'] * params_j['C6'])
                C12 = np.sqrt(params_i['C12'] * params_j['C12'])

                r6 = r**6
                r12 = r6 * r6

                e_lj += C12 / r12 - C6 / r6

                f_magnitude = (12 * C12 / (r12 * r * r) - 6 * C6 / (r6 * r * r))
                f_vec = f_magnitude * r_vec / r

                pi.force += f_vec
                pj.force -= f_vec

        # Bonds
        for bond in self.bonds:
            pi, pj = self.particles[bond.i], self.particles[bond.j]

            r_vec = pj.position - pi.position
            r_vec = self.apply_pbc(r_vec)
            r = np.linalg.norm(r_vec)

            if r < 0.001:
                continue

            dr = r - bond.r0
            e_bond += 0.5 * bond.k * dr * dr

            f_magnitude = -bond.k * dr
            f_vec = f_magnitude * r_vec / r

            pi.force += f_vec
            pj.force -= f_vec

        # Angles
        for angle in self.angles:
            pi = self.particles[angle.i]
            pj = self.particles[angle.j]
            pk = self.particles[angle.k]

            r_ij = pi.position - pj.position
            r_kj = pk.position - pj.position
            r_ij = self.apply_pbc(r_ij)
            r_kj = self.apply_pbc(r_kj)

            rij = np.linalg.norm(r_ij)
            rkj = np.linalg.norm(r_kj)

            if rij < 0.001 or rkj < 0.001:
                continue

            cos_theta = np.dot(r_ij, r_kj) / (rij * rkj)
            cos_theta = np.clip(cos_theta, -1.0, 1.0)
            theta = np.arccos(cos_theta)

            dtheta = theta - angle.theta0
            e_angle += 0.5 * angle.k_angle * dtheta * dtheta

            if abs(np.sin(theta)) > 0.001:
                f_factor = -angle.k_angle * dtheta / np.sin(theta)
                f_i = f_factor * (r_kj / (rij * rkj) - cos_theta * r_ij / (rij * rij))
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
        """Perform one integration step."""
        if self.integrator == 'velocity-verlet':
            self._velocity_verlet()
        elif self.integrator == 'leap-frog':
            self._leap_frog()

    def _velocity_verlet(self):
        """Velocity Verlet integration."""
        for p in self.particles:
            p.velocity += 0.5 * self.dt * (p.force / p.mass)

        for p in self.particles:
            p.position += self.dt * p.velocity
            if self.periodic:
                p.position = np.mod(p.position, self.box_size)

        energies = self.calculate_forces()

        for p in self.particles:
            p.velocity += 0.5 * self.dt * (p.force / p.mass)

        if self.target_temperature is not None and self.temperature_coupling is not None:
            self._berendsen_thermostat()

        return energies

    def _leap_frog(self):
        """Leap-Frog integration."""
        energies = self.calculate_forces()

        for p in self.particles:
            p.velocity += self.dt * (p.force / p.mass)

        for p in self.particles:
            p.position += self.dt * p.velocity
            if self.periodic:
                p.position = np.mod(p.position, self.box_size)

        if self.target_temperature is not None and self.temperature_coupling is not None:
            self._berendsen_thermostat()

        return energies

    def _berendsen_thermostat(self):
        """Berendsen thermostat."""
        current_temp = self.temperature()
        if current_temp < 0.001:
            return

        scale = np.sqrt(1 + self.dt / self.temperature_coupling *
                       (self.target_temperature / current_temp - 1))

        for p in self.particles:
            p.velocity *= scale

    def temperature(self) -> float:
        """Calculate instantaneous temperature (3D: 3 degrees of freedom)."""
        if len(self.particles) == 0:
            return 0.0

        ke = sum(p.kinetic_energy() for p in self.particles)
        k_B = 0.00831446261815324
        n_df = 3 * len(self.particles)

        return 2 * ke / (k_B * n_df) if n_df > 0 else 0.0

    def total_energy(self) -> Tuple[float, float, float]:
        """Calculate total energy."""
        ke = sum(p.kinetic_energy() for p in self.particles)
        energies = self.calculate_forces()
        pe = energies['potential']
        return ke, pe, ke + pe

    def run(self, n_steps: int, callback: Optional[Callable] = None):
        """Run simulation."""
        for _ in range(n_steps):
            energies = self.integrate_step()

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
        """Get all particle positions as numpy array (N, 3)."""
        return np.array([p.position for p in self.particles])

    def get_velocities(self) -> np.ndarray:
        """Get all particle velocities as numpy array (N, 3)."""
        return np.array([p.velocity for p in self.particles])

    def get_forces(self) -> np.ndarray:
        """Get all particle forces as numpy array (N, 3)."""
        return np.array([p.force for p in self.particles])

    def get_colors(self) -> List[str]:
        """Get particle colors for visualization."""
        return [p.color for p in self.particles]

    def set_temperature(self, temperature: float, tau: float = 0.1):
        """Enable temperature coupling."""
        self.target_temperature = temperature
        self.temperature_coupling = tau

    def initialize_velocities(self, temperature: float):
        """Initialize velocities from Maxwell-Boltzmann distribution (3D)."""
        k_B = 0.00831446261815324

        for p in self.particles:
            sigma = np.sqrt(k_B * temperature / p.mass)
            p.velocity = np.random.normal(0, sigma, size=3)

        # Remove center-of-mass motion
        total_momentum = sum(p.mass * p.velocity for p in self.particles)
        total_mass = sum(p.mass for p in self.particles)
        v_com = total_momentum / total_mass

        for p in self.particles:
            p.velocity -= v_com
