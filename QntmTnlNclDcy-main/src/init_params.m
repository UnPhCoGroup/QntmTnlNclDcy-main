function params = init_params(varargin)
%INIT_PARAMS Initialize simulation parameters for quantum tunneling simulation
%
%   This function sets up all physical constants, simulation parameters,
%   and numerical settings for the quantum tunneling simulation.
%
%   Inputs:
%       varargin - Optional parameter-value pairs to override defaults
%                  Supported parameters: 'barrier_type', 'barrier_height', 
%                  'barrier_width', 'well_width', 'particle_mass', 'energy'
%
%   Outputs:
%       params - Structure containing all simulation parameters
%                Fields: physics, numerical, potential, simulation
%
%   Example:
%       params = init_params('barrier_height', 50, 'energy', 30);
%
%   Author: UniPhi Collective
%   Date: July-September 2025

%% Physical Constants (SI units)
params.physics.hbar = 1.054571817e-34;  % Reduced Planck constant [J⋅s]
params.physics.e = 1.602176634e-19;     % Elementary charge [C]
params.physics.m_e = 9.1093837015e-31;  % Electron mass [kg]
params.physics.m_p = 1.67262192369e-27; % Proton mass [kg]
params.physics.m_alpha = 6.6446573357e-27; % Alpha particle mass [kg]
params.physics.k_B = 1.380649e-23;      % Boltzmann constant [J/K]

%% Default Simulation Parameters
params.simulation.barrier_type = 'square';  % 'square', 'double', or 'coulomb'
params.simulation.particle_type = 'alpha';  % 'alpha', 'electron', 'proton'

%% Potential Parameters
params.potential.barrier_height = 30;       % Barrier height [MeV]
params.potential.barrier_width = 2e-15;     % Barrier width [m] (2 fm)
params.potential.well_width = 5e-15;        % Well width [m] (5 fm)
params.potential.coulomb_radius = 7.5e-15;  % Nuclear radius [m] (7.5 fm - realistic for heavy nuclei)
params.potential.coulomb_charge = 2;        % Alpha particle charge [e]
params.potential.gap_width = 4e-15;         % Gap between double barriers [m] (4 fm)

%% Energy and Mass Parameters
params.energy.incident_energy = 4.2;        % Incident particle energy [MeV] (realistic for U-238 alpha decay)
params.energy.energy_range = [0.1, 50];     % Energy scan range [MeV]
params.energy.num_energies = 100;           % Number of energy points

% Set particle mass based on type
switch params.simulation.particle_type
    case 'alpha'
        params.physics.m = params.physics.m_alpha;
    case 'electron'
        params.physics.m = params.physics.m_e;
    case 'proton'
        params.physics.m = params.physics.m_p;
    otherwise
        params.physics.m = params.physics.m_alpha; % Default to alpha
end

%% Numerical Parameters
params.numerical.x_min = -20e-15;           % Left boundary [m]
params.numerical.x_max = 20e-15;            % Right boundary [m]
params.numerical.num_points = 2000;         % Number of grid points
params.numerical.tolerance = 1e-8;          % Eigenvalue convergence tolerance
params.numerical.max_iterations = 1000;     % Maximum iterations for solver

%% Derived Parameters
params.numerical.dx = (params.numerical.x_max - params.numerical.x_min) / ...
                      (params.numerical.num_points - 1);
params.numerical.x = linspace(params.numerical.x_min, ...
                              params.numerical.x_max, ...
                              params.numerical.num_points);

% Convert energy from MeV to Joules for calculations
params.energy.incident_energy_J = params.energy.incident_energy * ...
                                   params.physics.e * 1e6;
params.potential.barrier_height_J = params.potential.barrier_height * ...
                                     params.physics.e * 1e6;

%% Process Input Arguments
if nargin > 0
    for i = 1:2:length(varargin)
        param_name = varargin{i};
        param_value = varargin{i+1};
        
        % Handle nested parameter updates
        switch param_name
            case 'barrier_type'
                params.simulation.barrier_type = param_value;
            case 'gap_width'
                params.potential.gap_width = param_value;
            case 'barrier_height'
                params.potential.barrier_height = param_value;
                params.potential.barrier_height_J = param_value * ...
                                                    params.physics.e * 1e6;
            case 'barrier_width'
                params.potential.barrier_width = param_value;
            case 'well_width'
                params.potential.well_width = param_value;
            case 'particle_mass'
                params.physics.m = param_value;
            case 'energy'
                params.energy.incident_energy = param_value;
                params.energy.incident_energy_J = param_value * ...
                                                  params.physics.e * 1e6;
            case 'num_points'
                params.numerical.num_points = param_value;
                params.numerical.dx = (params.numerical.x_max - ...
                                      params.numerical.x_min) / ...
                                      (param_value - 1);
                params.numerical.x = linspace(params.numerical.x_min, ...
                                              params.numerical.x_max, ...
                                              param_value);
            otherwise
                warning('Unknown parameter: %s', param_name);
        end
    end
end

%% Validation
if params.potential.barrier_height <= 0
    error('Barrier height must be positive');
end

if params.energy.incident_energy <= 0
    error('Incident energy must be positive');
end

if params.numerical.num_points < 10
    error('Number of grid points must be at least 10');
end

%% Display Configuration Summary
fprintf('=== Quantum Tunneling Simulation Parameters ===\n');
fprintf('Particle: %s (mass = %.2e kg)\n', ...
        params.simulation.particle_type, params.physics.m);
fprintf('Barrier: %s, height = %.1f MeV, width = %.1f fm\n', ...
        params.simulation.barrier_type, ...
        params.potential.barrier_height, ...
        params.potential.barrier_width * 1e15);
fprintf('Energy: %.1f MeV\n', params.energy.incident_energy);
fprintf('Grid: %d points, dx = %.2e m\n', ...
        params.numerical.num_points, params.numerical.dx);
fprintf('===============================================\n');

end
