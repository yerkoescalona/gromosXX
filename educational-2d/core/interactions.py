"""
2D Interaction Potentials
==========================

Educational wrappers for gromos-rs interactions in 2D.
Demonstrates molecular force fields in a simplified 2D context.
"""

import numpy as np
from typing import Tuple, Optional


class LennardJones2D:
    """
    2D Lennard-Jones potential.

    The Lennard-Jones potential models van der Waals interactions:
    - Short-range repulsion (Pauli exclusion)
    - Long-range attraction (dispersion forces)

    GROMOS form: E(r) = C12/r^12 - C6/r^6

    Alternative form: E(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]

    Where:
    - epsilon: Depth of potential well
    - sigma: Distance where potential is zero
    - C12 = 4*epsilon*sigma^12
    - C6 = 4*epsilon*sigma^6
    """

    def __init__(self, epsilon: float = 1.0, sigma: float = 0.3):
        """
        Initialize Lennard-Jones parameters.

        Parameters
        ----------
        epsilon : float
            Well depth in kJ/mol
        sigma : float
            Zero-crossing distance in nm
        """
        self.epsilon = epsilon
        self.sigma = sigma

        # Convert to GROMOS C6/C12 form
        self.C6 = 4 * epsilon * sigma**6
        self.C12 = 4 * epsilon * sigma**12

    def energy(self, r: float) -> float:
        """
        Calculate LJ energy at distance r.

        Parameters
        ----------
        r : float
            Distance in nm

        Returns
        -------
        float
            Potential energy in kJ/mol
        """
        if r < 0.001:
            return 1e10  # Very large repulsion

        r6 = r**6
        r12 = r6 * r6

        return self.C12 / r12 - self.C6 / r6

    def force_magnitude(self, r: float) -> float:
        """
        Calculate LJ force magnitude at distance r.

        F(r) = -dE/dr = (12*C12/r^14 - 6*C6/r^8)

        Parameters
        ----------
        r : float
            Distance in nm

        Returns
        -------
        float
            Force magnitude (positive = repulsive)
        """
        if r < 0.001:
            return 1e10

        r6 = r**6
        r12 = r6 * r6
        r14 = r12 * r * r
        r8 = r6 * r * r

        return 12 * self.C12 / r14 - 6 * self.C6 / r8

    def minimum_energy(self) -> Tuple[float, float]:
        """
        Get minimum energy position and value.

        Returns
        -------
        Tuple[float, float]
            (r_min, E_min) where r_min = 2^(1/6) * sigma
        """
        r_min = 2**(1/6) * self.sigma
        e_min = -self.epsilon
        return r_min, e_min


class HarmonicBond2D:
    """
    2D harmonic bond potential.

    E(r) = 0.5 * k * (r - r0)^2

    Models covalent bonds with a spring-like potential.
    """

    def __init__(self, k: float = 1000.0, r0: float = 0.1):
        """
        Initialize harmonic bond parameters.

        Parameters
        ----------
        k : float
            Spring constant in kJ/(mol·nm²)
        r0 : float
            Equilibrium bond length in nm
        """
        self.k = k
        self.r0 = r0

    def energy(self, r: float) -> float:
        """
        Calculate bond energy.

        Parameters
        ----------
        r : float
            Current bond length in nm

        Returns
        -------
        float
            Bond energy in kJ/mol
        """
        dr = r - self.r0
        return 0.5 * self.k * dr * dr

    def force_magnitude(self, r: float) -> float:
        """
        Calculate bond force magnitude.

        F(r) = -k * (r - r0)

        Parameters
        ----------
        r : float
            Current bond length in nm

        Returns
        -------
        float
            Force magnitude (positive = repulsive)
        """
        return -self.k * (r - self.r0)


class AnglePotential2D:
    """
    2D angle bending potential.

    E(theta) = 0.5 * k * (theta - theta0)^2

    Models bond angle preferences (e.g., sp3 tetrahedral, sp2 planar).
    """

    def __init__(self, k: float = 100.0, theta0: float = np.pi):
        """
        Initialize angle parameters.

        Parameters
        ----------
        k : float
            Force constant in kJ/(mol·rad²)
        theta0 : float
            Equilibrium angle in radians
        """
        self.k = k
        self.theta0 = theta0

    def energy(self, theta: float) -> float:
        """
        Calculate angle bending energy.

        Parameters
        ----------
        theta : float
            Current angle in radians

        Returns
        -------
        float
            Angle energy in kJ/mol
        """
        dtheta = theta - self.theta0
        return 0.5 * self.k * dtheta * dtheta


class QuarticBond2D:
    """
    2D quartic bond potential (GROMOS default).

    E(r) = 0.25 * k * (r^2 - r0^2)^2

    More realistic than harmonic for large deviations.
    """

    def __init__(self, k: float = 1000.0, r0: float = 0.1):
        """
        Initialize quartic bond parameters.

        Parameters
        ----------
        k : float
            Force constant in kJ/(mol·nm⁴)
        r0 : float
            Equilibrium bond length in nm
        """
        self.k = k
        self.r0 = r0

    def energy(self, r: float) -> float:
        """Calculate quartic bond energy."""
        r2 = r * r
        r02 = self.r0 * self.r0
        diff = r2 - r02
        return 0.25 * self.k * diff * diff

    def force_magnitude(self, r: float) -> float:
        """
        Calculate quartic bond force.

        F(r) = -dE/dr = -k * r * (r^2 - r0^2)
        """
        if r < 0.001:
            return 0.0

        r2 = r * r
        r02 = self.r0 * self.r0
        return -self.k * r * (r2 - r02)


def combining_rules(params1: dict, params2: dict,
                   rule: str = 'lorentz-berthelot') -> dict:
    """
    Mixing rules for LJ parameters between different atom types.

    Parameters
    ----------
    params1, params2 : dict
        LJ parameters {'epsilon': ..., 'sigma': ...} or {'C6': ..., 'C12': ...}
    rule : str
        Combining rule: 'lorentz-berthelot' or 'geometric'

    Returns
    -------
    dict
        Combined parameters
    """
    if rule == 'lorentz-berthelot':
        # Lorentz-Berthelot: sigma_ij = (sigma_i + sigma_j)/2
        #                    epsilon_ij = sqrt(epsilon_i * epsilon_j)

        if 'sigma' in params1 and 'sigma' in params2:
            sigma = (params1['sigma'] + params2['sigma']) / 2
            epsilon = np.sqrt(params1['epsilon'] * params2['epsilon'])

            return {
                'sigma': sigma,
                'epsilon': epsilon,
                'C6': 4 * epsilon * sigma**6,
                'C12': 4 * epsilon * sigma**12
            }
        elif 'C6' in params1 and 'C6' in params2:
            # Direct C6/C12 combining
            C6 = np.sqrt(params1['C6'] * params2['C6'])
            C12 = np.sqrt(params1['C12'] * params2['C12'])

            # Back-calculate sigma and epsilon
            sigma = (C12 / C6)**(1/6) if C6 > 0 else 0
            epsilon = C6**2 / (4 * C12) if C12 > 0 else 0

            return {
                'C6': C6,
                'C12': C12,
                'sigma': sigma,
                'epsilon': epsilon
            }

    elif rule == 'geometric':
        # Geometric mean for both
        sigma = np.sqrt(params1['sigma'] * params2['sigma'])
        epsilon = np.sqrt(params1['epsilon'] * params2['epsilon'])

        return {
            'sigma': sigma,
            'epsilon': epsilon,
            'C6': 4 * epsilon * sigma**6,
            'C12': 4 * epsilon * sigma**12
        }

    else:
        raise ValueError(f"Unknown combining rule: {rule}")


def plot_potential(potential, r_range: Tuple[float, float] = (0.2, 1.0),
                  n_points: int = 100, title: str = "Potential Energy"):
    """
    Plot a potential energy curve.

    Parameters
    ----------
    potential : LennardJones2D, HarmonicBond2D, etc.
        The potential object
    r_range : Tuple[float, float]
        Range of distances to plot (r_min, r_max) in nm
    n_points : int
        Number of points to evaluate
    title : str
        Plot title
    """
    import matplotlib.pyplot as plt

    r = np.linspace(r_range[0], r_range[1], n_points)
    energies = [potential.energy(ri) for ri in r]
    forces = [potential.force_magnitude(ri) for ri in r]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Energy plot
    ax1.plot(r, energies, 'b-', linewidth=2)
    ax1.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax1.set_xlabel('Distance r (nm)')
    ax1.set_ylabel('Energy (kJ/mol)')
    ax1.set_title(f'{title} - Energy')
    ax1.grid(True, alpha=0.3)

    # Force plot
    ax2.plot(r, forces, 'r-', linewidth=2)
    ax2.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax2.set_xlabel('Distance r (nm)')
    ax2.set_ylabel('Force (kJ/mol/nm)')
    ax2.set_title(f'{title} - Force')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


# Example usage and educational demonstrations
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    print("=" * 60)
    print("Educational 2D Interaction Potentials")
    print("=" * 60)

    # 1. Lennard-Jones potential
    print("\n1. Lennard-Jones Potential")
    print("-" * 40)

    lj = LennardJones2D(epsilon=1.0, sigma=0.3)
    print(f"Parameters: epsilon = {lj.epsilon} kJ/mol, sigma = {lj.sigma} nm")
    print(f"GROMOS form: C6 = {lj.C6:.6f}, C12 = {lj.C12:.8f}")

    r_min, e_min = lj.minimum_energy()
    print(f"Minimum at r = {r_min:.4f} nm, E = {e_min:.4f} kJ/mol")

    plot_potential(lj, r_range=(0.25, 1.0), title="Lennard-Jones")

    # 2. Harmonic bond
    print("\n2. Harmonic Bond")
    print("-" * 40)

    bond = HarmonicBond2D(k=1000.0, r0=0.15)
    print(f"Parameters: k = {bond.k} kJ/(mol·nm²), r0 = {bond.r0} nm")

    plot_potential(bond, r_range=(0.1, 0.3), title="Harmonic Bond")

    # 3. Combining rules demonstration
    print("\n3. Combining Rules")
    print("-" * 40)

    params_C = {'epsilon': 0.5, 'sigma': 0.3}  # Carbon-like
    params_O = {'epsilon': 0.7, 'sigma': 0.25}  # Oxygen-like

    params_CO = combining_rules(params_C, params_O)
    print(f"C-C: epsilon = {params_C['epsilon']}, sigma = {params_C['sigma']}")
    print(f"O-O: epsilon = {params_O['epsilon']}, sigma = {params_O['sigma']}")
    print(f"C-O: epsilon = {params_CO['epsilon']:.4f}, sigma = {params_CO['sigma']:.4f}")

    plt.show()
