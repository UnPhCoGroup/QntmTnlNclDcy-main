# Quantum Tunneling In Nuclear Decay

This is a comprehensive MATLAB simulation of quantum tunneling in alpha decay, implementing the 1D time-independent Schrödinger equation using finite difference methods.

## Project Overview

This project simulates quantum tunneling phenomena, particularly focusing on alpha decay processes where alpha particles tunnel through Coulomb barriers. The simulation provides both theoretical understanding and practical visualization of quantum mechanical tunneling effects.

### Features

- **Triple Potential Support**: Square barrier, double barrier, and Coulomb barrier potentials
- **Finite Difference Solver**: Robust numerical solution of the Schrödinger equation
- **Transmission Analysis**: Complete transmission/reflection coefficient calculations
- **Comprehensive Visualization**: Static plots, animations, and parameter sweeps
- **Parameter Studies**: Systematic analysis of how parameters affect tunneling
- **Nuclear Physics Accuracy**: Realistic parameters for alpha decay (U-238, 4.2 MeV)
- **Smart Scaling**: Potential scaling for optimal visualization
- **Multiple Visualization Styles**: Traditional and scaled wavefunction plots

## Project Structure

```
src/
├── main.m                    # Main orchestrator function
├── init_params.m            # Parameter initialization and configuration
├── solve_schrodinger.m      # Core finite difference solver
├── potential_square.m       # Square barrier potential
├── potential_double_barrier.m # Double barrier potential (interference effects)
├── potential_coulomb.m      # Coulomb barrier potential (alpha decay)
├── compute_transmission.m   # Transmission coefficient calculations
├── visualize_wavefunction.m # Wavefunction visualization
├── visualize_single_wavefunction.m # James-style single wavefunction plots
├── animate_tunneling.m      # Animated tunneling visualization
├── animate_double_barrier_tunneling.m # Double barrier animation
├── sweep_parameters.m       # Parameter sweep analysis
├── data/                    # Data storage
│   ├── parameters/          # Configuration files
│   └── results/            # Simulation results
├── docs/                   # Documentation
└── figs/                   # Generated figures and animations
```

## Team Roles

- **Alex**: Core solver development (`solve_schrodinger.m`)
- **Michael**: Potential function implementation (`potential_square.m`, `potential_coulomb.m`)
- **Nate**: Transmission coefficient calculations (`compute_transmission.m`)
- **James**: Visualization and animation (`visualize_wavefunction.m`, `animate_tunneling.m`)

## Quick Start

### Basic Usage

```matlab
% Run simulation with default parameters
main();

% Run with custom parameters
main('barrier_height', 40, 'energy', 10, 'barrier_type', 'coulomb');
```

### Parameter Sweep Analysis

```matlab
% Sweep barrier height from 10 to 50 MeV
barrier_heights = 10:5:50;
results = sweep_parameters('barrier_height', barrier_heights);

% Sweep energy with fixed barrier
energies = 1:0.5:20;
results = sweep_parameters('energy', energies, 'barrier_height', 30);
```

### Animation

```matlab
% Create tunneling animation
params = init_params();
V = potential_square(params);
animate_tunneling(V, params, 'energy', 8, 'duration', 5);
```

## User Guide

### Core Physics

The simulation is based on the time-independent Schrödinger equation:

```
-ℏ²/(2m) d²ψ/dx² + V(x)ψ = Eψ
```

Where:
- `ψ(x)` is the wavefunction
- `V(x)` is the potential energy
- `E` is the energy eigenvalue
- `m` is the particle mass
- `ℏ` is the reduced Planck constant

The simulation uses finite difference methods to discretize this equation into a matrix eigenvalue problem:

```
Hψ = Eψ
```

Where `H` is the Hamiltonian matrix representing the kinetic and potential energy operators.

### Module Architecture

The simulation is organized into several specialized modules:

#### 1. **Parameter Management (`init_params.m`)**
- **Purpose**: Centralized configuration and parameter validation
- **Key Features**:
  - Physical constants (masses, charges, fundamental constants)
  - Numerical parameters (grid resolution, domain size)
  - Potential parameters (barrier heights, widths, types)
  - Energy and particle specifications
- **Usage**: `params = init_params('barrier_type', 'double', 'energy', 5);`

#### 2. **Potential Functions**
- **`potential_square.m`**: Single square barrier potential
- **`potential_double_barrier.m`**: Double barrier with configurable gap (interference effects)
- **`potential_coulomb.m`**: Realistic Coulomb barrier for alpha decay
- **Physics**: Each potential represents different physical scenarios:
  - Square barriers: Idealized tunneling studies
  - Double barriers: Quantum interference and resonance effects
  - Coulomb barriers: Nuclear alpha decay processes

#### 3. **Core Solver (`solve_schrodinger.m`)**
- **Method**: Finite difference discretization with sparse matrix techniques
- **Algorithm**: 
  1. Constructs tridiagonal Hamiltonian matrix
  2. Uses MATLAB's `eigs()` for eigenvalue computation
  3. Returns energy eigenvalues and corresponding wavefunctions
- **Validation**: Automatic normalization and boundary condition enforcement

#### 4. **Transmission Analysis (`compute_transmission.m`)**
- **Purpose**: Calculates transmission and reflection coefficients
- **Method**: Scattering theory approach with plane wave solutions
- **Output**: Energy-dependent transmission probabilities
- **Applications**: Tunneling probability studies, barrier effectiveness analysis

#### 5. **Visualization Suite**
- **`visualize_wavefunction.m`**: Comprehensive wavefunction analysis
- **`visualize_single_wavefunction.m`**: James-style scaled visualization
- **`animate_tunneling.m`**: Time-dependent wave packet evolution
- **`animate_double_barrier_tunneling.m`**: Specialized double barrier animations

#### 6. **Parameter Studies (`sweep_parameters.m`)**
- **Purpose**: Systematic parameter variation analysis
- **Capabilities**: Energy sweeps, barrier height studies, width optimization
- **Output**: Comprehensive data tables and visualization


#### Advanced Workflows

**Workflow 1: Double Barrier Interference Study**
```matlab
% Set up double barrier with specific gap
params = init_params('barrier_type', 'double', 'gap_width', 4e-15);

% Generate potential and solve
V = potential_double_barrier(params);
[eigenvalues, eigenfunctions] = solve_schrodinger(V, params);

% Analyze transmission
[transmission, reflection, energies] = compute_transmission(V, params);

% Visualize results
visualize_wavefunction(eigenvalues, eigenfunctions, V, params);
animate_double_barrier_tunneling(params, 'energy', 5, 'duration', 6);
```

**Workflow 2: Alpha Decay Analysis**
```matlab
% Realistic alpha decay parameters
params = init_params('barrier_type', 'coulomb', 'energy', 4.2, 'particle_type', 'alpha');

% Solve and analyze
V = potential_coulomb(params);
[eigenvalues, eigenfunctions] = solve_schrodinger(V, params);

% Compute Gamow factor and transmission
[transmission, reflection, energies] = compute_transmission(V, params);

% Create publication-quality plots
visualize_single_wavefunction(params.numerical.x, eigenfunctions, V, eigenvalues, 1);
```

**Workflow 3: Parameter Optimization**
```matlab
% Study effect of barrier height on tunneling
heights = 10:2:40;  % MeV
results = sweep_parameters('barrier_height', heights);

% Analyze results
figure;
semilogy(results.barrier_height, results.transmission_coefficient);
xlabel('Barrier Height (MeV)');
ylabel('Transmission Coefficient');
title('Tunneling vs Barrier Height');
```

### Understanding the Output

#### 1. **Energy Eigenvalues**
- **Units**: Joules (convert to MeV by dividing by 1.602e-19 × 1e6)
- **Interpretation**: Bound state energies for particles in the potential
- **Physical Meaning**: Quantized energy levels due to quantum confinement

#### 2. **Wavefunctions**
- **Normalization**: ∫|ψ(x)|²dx = 1
- **Symmetry**: Even/odd parity states for symmetric potentials
- **Penetration**: Exponential decay in classically forbidden regions

#### 3. **Transmission Coefficients**
- **Range**: 0 ≤ T ≤ 1
- **Physical Meaning**: Probability of particle transmission through barrier
- **Energy Dependence**: Resonant enhancement at specific energies

#### 4. **Visualization Elements**
- **Blue curves**: Probability density |ψ(x)|²
- **Red lines**: Potential energy V(x)
- **Green lines**: Energy eigenvalues E
- **Scaling**: James-style scaling for optimal visual comparison

### Parameter Tuning Guide

#### **Barrier Parameters**
- **`barrier_height`**: Controls tunneling difficulty (higher = less tunneling)
- **`barrier_width`**: Affects tunneling probability (wider = less tunneling)
- **`gap_width`**: For double barriers, controls interference effects

#### **Energy Parameters**
- **`incident_energy`**: Particle kinetic energy (MeV)
- **`energy_range`**: Sweep range for transmission studies
- **`num_energies`**: Resolution of energy sweeps

#### **Numerical Parameters**
- **`num_points`**: Grid resolution (higher = more accurate, slower)
- **`domain_size`**: Spatial extent of simulation
- **`boundary_condition`**: How wavefunction behaves at edges

#### **Particle Parameters**
- **`particle_type`**: 'alpha', 'electron', 'proton'
- **`mass`**: Particle mass (affects de Broglie wavelength)
- **`charge`**: For Coulomb interactions

### Troubleshooting

#### **Common Issues**

1. **"No bound states found"**
   - **Cause**: Potential too shallow or domain too small
   - **Solution**: Increase `barrier_height` or `domain_size`

2. **"Transmission coefficient > 1"**
   - **Cause**: Numerical instability or incorrect normalization
   - **Solution**: Increase `num_points` or check potential definition

3. **"Animation too slow"**
   - **Cause**: High resolution or long duration
   - **Solution**: Reduce `num_points`, `duration`, or `fps`

4. **"Unphysical results"**
   - **Cause**: Inconsistent units or unrealistic parameters
   - **Solution**: Check parameter values and units

#### **Performance Optimization**

- **Grid Resolution**: Start with `num_points = 1000`, increase as needed
- **Animation**: Use `fps = 15` for quick previews, `fps = 30` for final videos
- **Memory**: Close figures between runs with `close all;`

### Scientific Applications

#### **1. Nuclear Physics**
- Alpha decay lifetime calculations
- Gamow factor analysis
- Nuclear potential studies

#### **2. Solid State Physics**
- Semiconductor tunneling
- Quantum dot confinement
- Superconductor Josephson junctions

#### **3. Quantum Mechanics Education**
- Wave-particle duality demonstration
- Uncertainty principle visualization
- Quantum interference effects

#### **4. Research Applications**
- Parameter optimization for quantum devices
- Tunneling probability calculations
- Resonance condition analysis

### Best Practices

#### **1. Parameter Selection**
- Start with default parameters
- Make incremental changes
- Validate results against known cases
- Use appropriate units (nuclear vs SI)

#### **2. Visualization**
- Use James-style scaling for publication figures
- Include error bars for parameter studies
- Save figures in multiple formats (.png, .fig)
- Add comprehensive legends and labels

#### **3. Analysis**
- Always check normalization
- Verify energy conservation
- Compare with analytical solutions when possible
- Document parameter choices

#### **4. Code Organization**
- Use descriptive variable names
- Comment complex calculations
- Save intermediate results
- Version control parameter files

### Example Research Projects

#### **Project 1: Double Barrier Resonance Study**
```matlab
% Investigate resonance conditions in double barriers
gap_widths = 2e-15:0.5e-15:8e-15;
resonance_energies = zeros(size(gap_widths));

for i = 1:length(gap_widths)
    params = init_params('barrier_type', 'double', 'gap_width', gap_widths(i));
    V = potential_double_barrier(params);
    [eigenvalues, ~] = solve_schrodinger(V, params);
    resonance_energies(i) = eigenvalues(1) / (1.602e-19 * 1e6);  % Convert to MeV
end

plot(gap_widths * 1e15, resonance_energies, 'o-');
xlabel('Gap Width (fm)');
ylabel('Resonance Energy (MeV)');
title('Double Barrier Resonance vs Gap Width');
```

#### **Project 2: Alpha Decay Lifetime Calculation**
```matlab
% Calculate alpha decay lifetime for different nuclei
nuclei = {'U-238', 'Th-232', 'Ra-226'};
charges = [92, 90, 88];
energies = [4.2, 4.0, 4.8];  % MeV

lifetimes = zeros(size(nuclei));

for i = 1:length(nuclei)
    params = init_params('barrier_type', 'coulomb', 'energy', energies(i));
    params.potential.coulomb_charge = charges(i);
    V = potential_coulomb(params);
    [transmission, ~, ~] = compute_transmission(V, params);
    
    % Calculate lifetime (simplified)
    lifetimes(i) = 1 / transmission;
end

semilogy(1:length(nuclei), lifetimes, 'o-');
set(gca, 'XTick', 1:length(nuclei), 'XTickLabel', nuclei);
ylabel('Relative Lifetime');
title('Alpha Decay Lifetime Comparison');
```

## Configuration

### Available Parameters

| Parameter | Description | Default | Units |
|-----------|-------------|---------|-------|
| `barrier_type` | Type of potential barrier | `'square'` | - |
| `barrier_height` | Height of potential barrier | 30 | MeV |
| `barrier_width` | Width of square barrier | 2e-15 | m |
| `well_width` | Width of potential well | 5e-15 | m |
| `gap_width` | Gap between double barriers | 4e-15 | m |
| `energy` | Incident particle energy | 5 | MeV |
| `particle_type` | Type of particle | `'alpha'` | - |
| `num_points` | Grid resolution | 2000 | - |

### Physical Constants

The simulation uses accurate physical constants:
- Reduced Planck constant: ℏ = 1.054571817×10⁻³⁴ J⋅s
- Elementary charge: e = 1.602176634×10⁻¹⁹ C
- Alpha particle mass: m_α = 6.6446573357×10⁻²⁷ kg

## Scientific Background

### Quantum Tunneling

Quantum tunneling is a fundamental quantum mechanical phenomenon where particles can pass through potential energy barriers that would be classically insurmountable. This occurs due to the wave-like nature of particles described by the Schrödinger equation.

### Alpha Decay

Alpha decay is a nuclear decay process where an alpha particle (helium nucleus) escapes from a parent nucleus by tunneling through the Coulomb barrier. The probability of decay is determined by:

1. **Gamow Factor**: Quantifies the tunneling probability through the Coulomb barrier
2. **Nuclear Potential**: Describes the binding energy within the nucleus
3. **Energy Levels**: Determines the available decay channels

### Mathematical Framework

The time-independent Schrödinger equation:

```
-ℏ²/(2m) * d²ψ/dx² + V(x)ψ(x) = Eψ(x)
```

Where:
- ℏ is the reduced Planck constant
- m is the particle mass
- V(x) is the potential energy
- ψ(x) is the wavefunction
- E is the energy eigenvalue

## Output Files

### Figures
- `wavefunction_analysis.png`: Combined wavefunction analysis
- `detailed_wavefunctions.png`: Detailed wavefunction plots
- `energy_levels.png`: Energy level diagram
- `transmission_coefficients.png`: Transmission vs energy plots
- `parameter_sweep_analysis.png`: Parameter sweep results

### Data Files
- `wavefunction_profiles.mat`: Eigenvalues and eigenfunctions
- `transmission_vs_energy.mat`: Transmission coefficient data
- `parameter_sweep.mat`: Parameter sweep results

### Animations
- `tunneling_animation.mp4`: Animated tunneling process

## Requirements

- MATLAB R2023a or later
- Image Processing Toolbox (for video creation)
- Signal Processing Toolbox (recommended)

## Examples

### Example 1: Basic Square Barrier

```matlab
% Set up square barrier simulation
params = init_params('barrier_type', 'square', 'barrier_height', 25, 'energy', 8);
V = potential_square(params);
[eigenvalues, eigenfunctions] = solve_schrodinger(V, params);
visualize_wavefunction(eigenvalues, eigenfunctions, V, params);
```

### Example 1b: James-style Single Wavefunction

```matlab
% Create James-style scaled visualization
params = init_params();
V = potential_square(params);
[eigenvalues, eigenfunctions] = solve_schrodinger(V, params);
visualize_single_wavefunction(params.numerical.x, eigenfunctions, V, eigenvalues, 1);
```

### Example 2: Double Barrier (Reference Style)

```matlab
% Set up double barrier simulation (like the reference)
params = init_params('barrier_type', 'double', 'gap_width', 4e-15, 'energy', 5);
V = potential_double_barrier(params);
[eigenvalues, eigenfunctions] = solve_schrodinger(V, params);
visualize_wavefunction(eigenvalues, eigenfunctions, V, params);

% Create double barrier animation
animate_double_barrier_tunneling(params, 'energy', 5, 'duration', 6);
```

### Example 3: Coulomb Barrier (Alpha Decay)

```matlab
% Set up alpha decay simulation
params = init_params('barrier_type', 'coulomb', 'energy', 5);
V = potential_coulomb(params);
[transmission, reflection, energies] = compute_transmission(V, params);
plot(energies, transmission);
```

### Example 4: Parameter Study

```matlab
% Study effect of barrier height on tunneling
heights = 10:2:40;
results = sweep_parameters('barrier_height', heights);
plot(results.param_values, results.transmission_max);
```

## Validation and Testing

The simulation includes comprehensive validation:
- Energy conservation checks
- Probability normalization verification
- Orthogonality testing for eigenfunctions
- Physical bounds checking (transmission coefficients 0-1)

## Future Enhancements

- 3D potential support
- Time-dependent simulations
- Multi-particle systems
- Experimental data comparison
- Machine learning parameter optimization

## References

1. Griffiths, D.J. "Introduction to Quantum Mechanics" (3rd Edition)
2. Newman, M.E.J. "Computational Physics" 
3. Krane, K.S. "Introductory Nuclear Physics"
4. Shankar, R. "Principles of Quantum Mechanics"

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

This is a collaborative project. Please follow the established coding standards:
- Use consistent function signatures
- Include comprehensive documentation
- Validate all inputs and outputs
- Follow MATLAB best practices for vectorization

## Contact

For questions or contributions, please contact the development team.

---

*Developed by UniPhi Collective - July-September 2025*
