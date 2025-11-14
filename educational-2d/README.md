# Educational 2D Molecular Dynamics Simulations

**Interactive 2D MD simulations using GROMOS-RS refined interactions for teaching and learning molecular dynamics**

---

## Overview

This educational package provides simplified 2D molecular dynamics simulations that demonstrate fundamental concepts in molecular simulation using the refined force field interactions from **gromos-rs**. By working in 2D rather than 3D, students and researchers can more easily visualize and understand:

- Molecular forces and potentials
- Energy conservation and conversion
- Temperature and thermodynamics
- Structural transitions (gas → liquid → solid)
- Polymer physics
- Protein folding basics

## Features

✨ **Core Simulation Engine**
- 2D projection of GROMOS force field interactions
- Lennard-Jones (van der Waals) interactions
- Harmonic bonds and angle potentials
- Velocity Verlet and Leap-Frog integrators
- Berendsen thermostat for temperature control
- Periodic boundary conditions

🎨 **Real-Time Visualization**
- Matplotlib-based rendering
- Particle trajectories (trails)
- Energy plots
- Force/velocity vectors
- Customizable colors and styles

📚 **Educational Examples**
1. Gas in a box (ideal gas behavior)
2. Liquid droplet formation (phase separation)
3. Crystal lattice (solid-liquid-gas transitions)
4. Polymer chain dynamics (conformational sampling)
5. Protein folding (hydrophobic collapse)

📓 **Interactive Notebooks**
- Jupyter notebooks with step-by-step explanations
- Interactive parameter exploration
- Hands-on exercises

---

## Installation

### Prerequisites

```bash
# Python 3.8+
# Required packages:
pip install numpy matplotlib jupyter
```

### Optional: GROMOS-RS Python Bindings

For full integration with gromos-rs (recommended but not required):

```bash
cd ../py-gromos
pip install -e .
```

The educational simulations use standalone implementations but can leverage gromos-rs Python bindings when available.

---

## Quick Start

### Running Examples

```bash
cd examples/

# Example 1: Gas in a Box
python 01_gas_in_box.py

# Example 2: Liquid Droplet
python 02_liquid_droplet.py

# Example 3: Crystal Lattice
python 03_crystal_lattice.py

# Example 4: Polymer Chain
python 04_polymer_chain.py

# Example 5: Protein Folding
python 05_protein_folding_2d.py
```

### Interactive Notebooks

```bash
cd notebooks/
jupyter notebook 01_interactive_lennard_jones.ipynb
```

### Creating Your Own Simulation

```python
from core import Simulation2D, Renderer2D

# Create simulation
sim = Simulation2D(
    box_size=(5.0, 5.0),  # nm
    periodic=True,
    integrator='velocity-verlet',
    dt=0.002  # ps
)

# Add particles
for i in range(50):
    sim.add_particle(
        position=[np.random.uniform(0, 5), np.random.uniform(0, 5)],
        mass=1.0,
        color='blue'
    )

# Initialize velocities at 300 K
sim.initialize_velocities(temperature=300.0)
sim.set_temperature(300.0, tau=0.1)

# Create renderer
renderer = Renderer2D(box_size=(5.0, 5.0))

# Run simulation
for step in range(1000):
    sim.integrate_step()
    if step % 10 == 0:
        renderer.render(sim)

renderer.show()
```

---

## Educational Examples

### 1. Gas in a Box

**Demonstrates:** Ideal gas behavior, Maxwell-Boltzmann distribution, energy equipartition

```bash
python examples/01_gas_in_box.py
```

**Learning Objectives:**
- Understand kinetic theory of gases
- Observe temperature as emergent property
- See energy conservation in isolated system
- Verify equipartition theorem (⟨KE⟩ = kᵦT per degree of freedom)

**Key Output:**
- Real-time particle motion
- Energy vs time plots
- Temperature evolution
- Theoretical vs observed kinetic energy

---

### 2. Liquid Droplet Formation

**Demonstrates:** Phase separation, surface tension, cohesive forces

```bash
python examples/02_liquid_droplet.py
```

**Learning Objectives:**
- Observe liquid droplet coalescence
- Understand role of attractive interactions
- See surface tension in action
- Measure radius of gyration

**Key Output:**
- Two clusters merging into one
- Particle trajectories showing collective motion
- Radial distribution function
- Energy minimization

---

### 3. Crystal Lattice

**Demonstrates:** Solid-liquid transition, melting, phonons, order-disorder

```bash
python examples/03_crystal_lattice.py
```

**Learning Objectives:**
- Perfect hexagonal lattice formation
- Thermal vibrations (phonons) in solid
- Melting as temperature increases
- Order parameter evolution

**Key Output:**
- Three phases: solid (ordered) → melting → liquid (disordered)
- Mean square displacement from lattice
- Energy changes during phase transition
- Visual demonstration of melting

---

### 4. Polymer Chain Dynamics

**Demonstrates:** Bonded interactions, conformational sampling, polymer physics

```bash
python examples/04_polymer_chain.py
```

**Learning Objectives:**
- Covalent bonds as spring potentials
- Bond angle preferences
- Radius of gyration (Rg)
- End-to-end distance (Ree)
- Comparison with random walk theory

**Key Output:**
- Flexible chain exploring conformational space
- Rg and Ree time series
- Conformational statistics
- Comparison with theoretical predictions

---

### 5. Protein Folding (2D HP Model)

**Demonstrates:** Hydrophobic effect, protein folding, native state

```bash
python examples/05_protein_folding_2d.py
```

**Learning Objectives:**
- Hydrophobic-Polar (HP) model
- Hydrophobic collapse
- Native contact formation
- Folding funnel concept

**Key Output:**
- Extended → collapsed transition
- Hydrophobic residues cluster in core
- Polar residues on surface
- Contact map and folding analysis

---

## Core Modules

### `core/simulation2d.py`

Main simulation engine with:
- `Simulation2D`: 2D MD simulation class
- `Particle`: Particle data structure
- `Bond`: Harmonic bond definition
- `Angle`: Angle potential definition

**Key Methods:**
- `add_particle()`: Add particle to system
- `add_bond()`: Create covalent bond
- `add_angle()`: Add angle potential
- `integrate_step()`: Perform one MD step
- `calculate_forces()`: Compute all forces
- `temperature()`: Calculate instantaneous temperature
- `run()`: Execute simulation

---

### `core/renderer.py`

Visualization and rendering:
- `Renderer2D`: Real-time matplotlib visualization
- `plot_energy_summary()`: Energy analysis plots

**Features:**
- Particle rendering with custom colors
- Bond visualization
- Force/velocity vectors
- Trajectory trails
- Multi-panel energy plots

---

### `core/interactions.py`

Interaction potential wrappers:
- `LennardJones2D`: LJ potential (ε, σ parameterization)
- `HarmonicBond2D`: Harmonic spring (k, r₀)
- `AnglePotential2D`: Angle bending (k, θ₀)
- `QuarticBond2D`: Quartic bond (GROMOS style)
- `combining_rules()`: LJ parameter mixing
- `plot_potential()`: Visualize potential curves

---

## Physical Background

### Lennard-Jones Potential

The fundamental nonbonded interaction:

```
E(r) = 4ε[(σ/r)¹² - (σ/r)⁶]
     = C₁₂/r¹² - C₆/r⁶  (GROMOS form)
```

- **r⁻¹² term**: Short-range repulsion (Pauli exclusion)
- **r⁻⁶ term**: Long-range attraction (van der Waals)
- **ε**: Well depth (interaction strength)
- **σ**: Zero-crossing distance (particle size)
- **r_min = 2^(1/6)σ**: Minimum energy distance

### Harmonic Bond

Covalent bonds modeled as springs:

```
E(r) = ½k(r - r₀)²
```

- **k**: Spring constant (bond stiffness)
- **r₀**: Equilibrium bond length
- **F = -k(r - r₀)**: Restoring force

### Angle Potential

Bond angle preferences:

```
E(θ) = ½k_angle(θ - θ₀)²
```

- **θ**: Current angle between three atoms
- **θ₀**: Equilibrium angle
- **k_angle**: Angle stiffness

### Integration Algorithms

**Velocity Verlet** (default, more accurate):
```
v(t+½dt) = v(t) + ½a(t)dt
r(t+dt)  = r(t) + v(t+½dt)dt
a(t+dt)  = F(t+dt)/m
v(t+dt)  = v(t+½dt) + ½a(t+dt)dt
```

**Leap-Frog** (GROMOS standard, symplectic):
```
v(t+½dt) = v(t-½dt) + a(t)dt
r(t+dt)  = r(t) + v(t+½dt)dt
```

### Temperature Control

**Berendsen Thermostat** (velocity rescaling):
```
λ = √[1 + (dt/τ)(T₀/T - 1)]
v → λv
```

- **τ**: Coupling time constant
- **T₀**: Target temperature
- **T**: Current temperature

---

## Directory Structure

```
educational-2d/
├── README.md                    # This file
├── core/                        # Core simulation engine
│   ├── __init__.py
│   ├── simulation2d.py         # Main MD simulator
│   ├── renderer.py             # Visualization
│   └── interactions.py         # Force field potentials
├── examples/                    # Educational examples
│   ├── 01_gas_in_box.py
│   ├── 02_liquid_droplet.py
│   ├── 03_crystal_lattice.py
│   ├── 04_polymer_chain.py
│   └── 05_protein_folding_2d.py
├── notebooks/                   # Jupyter notebooks
│   └── 01_interactive_lennard_jones.ipynb
└── assets/                      # Output images/animations
    └── (generated during simulations)
```

---

## Tips for Educators

### Classroom Use

1. **Start with notebooks**: Interactive exploration builds intuition
2. **Run examples live**: Show real-time dynamics during lectures
3. **Modify parameters**: Let students experiment with different settings
4. **Compare with theory**: Use analytical predictions to verify results
5. **Assign projects**: Students can create their own simulation scenarios

### Suggested Progression

1. **Week 1**: Lennard-Jones potential (notebook + gas example)
2. **Week 2**: Bonded interactions (polymer example)
3. **Week 3**: Phase transitions (liquid/crystal examples)
4. **Week 4**: Protein folding (HP model example)
5. **Final Project**: Students design custom simulation

### Discussion Questions

- Why does increasing temperature melt a crystal?
- How does bond stiffness affect polymer conformations?
- What drives protein folding in the HP model?
- How do we know energy is conserved?
- What happens if we change the LJ parameters?

---

## Connections to GROMOS-RS

This educational package demonstrates the same physics used in production MD simulations:

| Educational 2D | GROMOS-RS 3D |
|----------------|--------------|
| `LennardJones2D` | `gromos-rs/src/interaction/nonbonded.rs` |
| `HarmonicBond2D` | `gromos-rs/src/interaction/bonded.rs` |
| `AnglePotential2D` | `gromos-rs/src/interaction/bonded.rs` |
| `Velocity Verlet` | `gromos-rs/src/integrator.rs` |
| `Berendsen thermostat` | `gromos-rs/src/algorithm/thermostats.rs` |

**Key Differences:**
- 2D vs 3D (simpler visualization)
- Simplified electrostatics (no PME/Ewald)
- Educational focus (clarity over performance)
- Standalone code (no dependencies on gromos-rs)

---

## Performance Notes

These simulations are designed for **education**, not production:

- **N < 100 particles**: Interactive real-time visualization
- **N = 100-500**: Moderate speed, good for demonstrations
- **N > 500**: May be slow due to matplotlib rendering

For large-scale simulations, use the full **gromos-rs** engine.

---

## Troubleshooting

### Simulation is unstable (particles flying apart)

- **Reduce timestep**: Try `dt=0.0005` instead of `dt=0.002`
- **Check initial positions**: Avoid particle overlaps
- **Increase temperature coupling**: Lower `tau` in `set_temperature()`

### Visualization is slow

- **Reduce render frequency**: `if step % 20 == 0: renderer.render(sim)`
- **Fewer particles**: Start with N=20-50
- **Disable trails**: `show_trails=False`

### Energy is not conserved

- **This is normal with thermostat**: Berendsen adds/removes energy
- **Without thermostat**: Check timestep is small enough
- **Check force calculations**: Verify LJ parameters are reasonable

---

## Contributing

We welcome contributions! Ideas for new examples:

- Diffusion coefficient calculation
- Radial distribution functions
- Hydrogen bonding (2D representation)
- Membrane formation (amphiphilic molecules)
- Enzyme-substrate binding
- Allostery and conformational change

Please follow the existing code style and add educational documentation.

---

## References

### Molecular Dynamics

- **Frenkel & Smit**: *Understanding Molecular Simulation* (2002)
- **Allen & Tildesley**: *Computer Simulation of Liquids* (2017)
- **Leach**: *Molecular Modelling: Principles and Applications* (2001)

### GROMOS Force Field

- **van Gunsteren et al.**: *Biomolecular Simulation: The GROMOS96 Manual and User Guide* (1996)
- **Oostenbrink et al.**: *J. Comput. Chem.* **25**, 1656 (2004)

### Educational Protein Models

- **Dill**: *Biochemistry* **24**, 1501 (1985) - HP model
- **Clementi et al.**: *J. Mol. Biol.* **298**, 937 (2000) - Structure-based models

---

## License

This educational package is part of the GROMOS-XX project and is distributed under the GPL-2.0 license.

---

## Acknowledgments

Built with:
- **gromos-rs**: High-performance Rust MD engine
- **matplotlib**: Python visualization
- **NumPy**: Numerical computing

Inspired by educational MD codes:
- **LAMMPS tutorials**
- **MDAnalysis examples**
- **PyMOL educational resources**

---

## Contact

For questions, suggestions, or educational collaborations:
- Open an issue on the GROMOS-XX GitHub repository
- Contribute examples or improvements via pull request

**Happy Simulating! 🎓⚛️**
