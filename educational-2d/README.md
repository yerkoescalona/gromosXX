# Educational 2D/3D Molecular Dynamics Simulations

**Interactive 2D/3D MD simulations with movie rendering using GROMOS-RS refined interactions for teaching and learning molecular dynamics**

---

## Overview

This educational package provides simplified 2D and 3D molecular dynamics simulations that demonstrate fundamental concepts in molecular simulation using the refined force field interactions from **gromos-rs**. By working in 2D and 3D with easy visualization and movie creation, students and researchers can:

- **Visualize** molecular forces and potentials in real-time
- **Understand** energy conservation and conversion
- **Explore** temperature and thermodynamics
- **Observe** structural transitions (gas → liquid → solid)
- **Study** polymer physics and protein folding
- **Create** movies and animations for presentations
- **Compare** 2D vs 3D behavior side-by-side

## Features

✨ **Core Simulation Engine (2D & 3D)**
- Full GROMOS force field interactions
- Lennard-Jones (van der Waals) interactions
- Harmonic bonds and angle potentials
- Velocity Verlet and Leap-Frog integrators
- Berendsen thermostat for temperature control
- Periodic boundary conditions
- **NEW: 3D simulation support with same physics**

🎨 **Real-Time Visualization**
- Matplotlib-based 2D/3D rendering
- Particle trajectories (trails)
- Energy plots
- Force/velocity vectors
- Customizable colors and styles
- **NEW: Interactive 3D rotation and zoom**
- **NEW: Rotating camera views**

🎬 **Movie Rendering**
- **MP4 video export** (ffmpeg)
- **Animated GIF creation** (Pillow)
- **Frame sequence export**
- **Rotating 3D views**
- **Side-by-side comparisons**
- Easy one-line movie creation

📚 **Educational Examples (2D)**
1. Gas in a box (ideal gas behavior)
2. Liquid droplet formation (phase separation)
3. Crystal lattice (solid-liquid-gas transitions)
4. Polymer chain dynamics (conformational sampling)
5. Protein folding (hydrophobic collapse)

🔮 **Educational Examples (3D)**
6. Gas in a cube (3D ideal gas with rotation)
7. Water-like molecule (molecular vibrations in 3D)

📓 **Interactive Notebooks**
- Jupyter notebooks with step-by-step explanations
- Interactive parameter exploration
- Hands-on exercises
- Lennard-Jones potential deep-dive

---

## Installation

### Prerequisites

```bash
# Python 3.8+
# Required packages:
pip install numpy matplotlib jupyter

# For movie creation (optional but recommended):
pip install Pillow  # For GIF animations

# For MP4 videos (requires ffmpeg):
conda install -c conda-forge ffmpeg
# or on Ubuntu/Debian:
# sudo apt-get install ffmpeg
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

# 2D Examples
python 01_gas_in_box.py          # Gas in a box
python 02_liquid_droplet.py      # Liquid droplet formation
python 03_crystal_lattice.py     # Crystal lattice
python 04_polymer_chain.py       # Polymer chain dynamics
python 05_protein_folding_2d.py  # Protein folding (HP model)

# 3D Examples
python 06_3d_gas_cube.py         # 3D gas with rotation
python 07_3d_molecule.py         # Water-like molecule vibrations

# Movie Creation
python movie_demo_2d.py          # Create MP4s and GIFs
```

### Interactive Notebooks

```bash
cd notebooks/
jupyter notebook 01_interactive_lennard_jones.ipynb
```

### Creating Movies (One-Liner!)

```python
from core import Simulation2D, Renderer2D, make_quick_movie

# Setup simulation
sim = Simulation2D(box_size=(5.0, 5.0))
renderer = Renderer2D(box_size=(5.0, 5.0))

# ... add particles, initialize ...

# Create movie in one line!
make_quick_movie(sim, renderer, n_steps=1000, filename='my_simulation.mp4')
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

## NEW: Movie Creation 🎬

### Quick Movie Creation

Create publication-quality movies with a single function call:

```python
from core import make_quick_movie, make_gif_animation, make_high_quality_movie

# MP4 video (requires ffmpeg)
make_quick_movie(sim, renderer, n_steps=1000, filename='simulation.mp4')

# Animated GIF (requires Pillow)
make_gif_animation(sim, renderer, n_steps=500, filename='animation.gif')

# High-quality video (60 FPS, high DPI)
make_high_quality_movie(sim, renderer, n_steps=500, filename='hq.mp4')
```

### Advanced Movie Options

```python
from core.movie_utils import record_simulation_movie

record_simulation_movie(
    sim, renderer,
    n_steps=2000,
    output_file='custom.mp4',
    render_interval=2,      # Render every 2 steps (faster)
    fps=30,                  # Frames per second
    dpi=150,                 # Resolution
    show_bonds=True,
    progress_interval=100
)
```

### Rotating 3D Movies

```python
from core.movie_utils import make_rotating_3d_movie

# Create 360° rotating view of 3D system
make_rotating_3d_movie(
    sim3d, renderer3d,
    n_frames=360,
    output_file='rotation.mp4',
    elevation=30,
    fps=30
)
```

### Movie Formats

| Format | Use Case | Requirements |
|--------|----------|--------------|
| **MP4** | Presentations, video players | ffmpeg |
| **GIF** | Websites, documentation | Pillow |
| **Frames** | Custom post-processing | None |

---

## NEW: 3D Simulations 🔮

### 3D Gas Simulation

```python
from core import Simulation3D, Renderer3D

# Create 3D simulation
sim = Simulation3D(
    box_size=(5.0, 5.0, 5.0),  # Cubic box
    periodic=True,
    integrator='velocity-verlet',
    dt=0.002
)

# Add particles in 3D
for i in range(50):
    sim.add_particle(position=[x, y, z], color='blue')

sim.initialize_velocities(300.0)  # 300 K

# Create 3D renderer
renderer = Renderer3D(
    box_size=(5.0, 5.0, 5.0),
    azimuth=45,     # Viewing angle
    elevation=30
)

# Run with visualization
for step in range(1000):
    sim.integrate_step()

    if step % 10 == 0:
        # Optionally rotate view during simulation
        azimuth = 45 + (360 * step / 1000)
        renderer.render(sim, azimuth=azimuth)
```

### Interactive 3D Rendering

```python
from core.renderer3d import InteractiveRenderer3D

# Create interactive renderer
renderer = InteractiveRenderer3D(box_size=(5.0, 5.0, 5.0))

# Run simulation
for step in range(1000):
    sim.integrate_step()
    if step % 5 == 0:
        renderer.render(sim)

# User can:
# - Click and drag to rotate
# - Scroll to zoom
# - Fully interactive matplotlib 3D view
```

### 3D vs 2D Comparison

| Feature | 2D | 3D |
|---------|----|----|
| Degrees of freedom | 2 per particle | 3 per particle |
| Equipartition | ⟨KE⟩ = k_B T | ⟨KE⟩ = 3/2 k_B T |
| Visualization | Easy to see all particles | Requires rotation |
| Performance | Faster | Slightly slower |
| Realism | Educational | More realistic |

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

### `core/simulation3d.py`

3D simulation engine (same physics, 3D space):
- `Simulation3D`: 3D MD simulation class
- `Particle3D`: 3D particle structure
- `Bond3D`: 3D bonds
- `Angle3D`: 3D angles

**Same interface as 2D:**
- All methods identical to `Simulation2D`
- Just use 3D positions `[x, y, z]` instead of `[x, y]`
- Temperature calculation accounts for 3 DOF

---

### `core/renderer3d.py`

3D visualization:
- `Renderer3D`: Standard 3D rendering with rotating view
- `InteractiveRenderer3D`: Mouse-controlled rotation and zoom

**Features:**
- Matplotlib 3D projection
- Rotatable camera (azimuth, elevation)
- Interactive mouse controls
- Box edge visualization
- Energy plots (2D subplot)

---

### `core/movie_utils.py`

Movie creation utilities:
- `MovieRecorder`: Frame capture and export class
- `record_simulation_movie()`: General movie recording
- `make_quick_movie()`: Fast MP4 creation
- `make_high_quality_movie()`: High FPS/DPI video
- `make_gif_animation()`: Animated GIF export
- `make_rotating_3d_movie()`: 360° rotation movie

**Supported Formats:**
- MP4 (via ffmpeg)
- GIF (via Pillow)
- PNG frame sequences

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
