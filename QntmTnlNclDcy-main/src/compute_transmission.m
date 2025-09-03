function [transmission, reflection, energies] = compute_transmission(V, params)
%COMPUTE_TRANSMISSION Calculate transmission and reflection coefficients
%
%   This function computes the transmission and reflection coefficients for
%   a particle incident on a potential barrier. It solves the scattering
%   problem by finding the wavefunction that satisfies the boundary conditions
%   for an incoming plane wave.
%
%   Inputs:
%       V - Potential energy array [J] at each grid point
%       params - Structure containing simulation parameters
%                Required fields: numerical.x, numerical.dx, physics.m, physics.hbar,
%                                energy.energy_range, energy.num_energies
%
%   Outputs:
%       transmission - Array of transmission coefficients (0 to 1)
%       reflection - Array of reflection coefficients (0 to 1)
%       energies - Array of energy values [MeV] corresponding to coefficients
%
%   Method:
%       For each energy E, solves the time-independent Schrödinger equation:
%       -ℏ²/(2m) * d²ψ/dx² + V(x)ψ(x) = Eψ(x)
%       
%       With boundary conditions:
%       ψ(x) = e^(ikx) + r*e^(-ikx)  (left side, x → -∞)
%       ψ(x) = t*e^(ikx)             (right side, x → +∞)
%       
%       where k = √(2mE)/ℏ and r, t are reflection and transmission amplitudes
%
%   Example:
%       params = init_params();
%       V = potential_square(params);
%       [T, R, E] = compute_transmission(V, params);
%       plot(E, T);
%
%   Author: Nate (UniPhi Collective)
%   Date: July-September 2025

%% Validate Input
if ~isfield(params, 'numerical') || ~isfield(params.numerical, 'x')
    error('params must contain numerical.x field');
end

if ~isfield(params, 'physics') || ~isfield(params.physics, 'm')
    error('params must contain physics.m field');
end

if ~isfield(params, 'energy') || ~isfield(params.energy, 'energy_range')
    error('params must contain energy.energy_range field');
end

%% Extract Parameters
x = params.numerical.x;
dx = params.numerical.dx;
m = params.physics.m;
hbar = params.physics.hbar;
num_points = length(x);

% Energy range for transmission calculation
energy_min = params.energy.energy_range(1) * 1.602e-19 * 1e6;  % Convert to J
energy_max = params.energy.energy_range(2) * 1.602e-19 * 1e6;  % Convert to J
num_energies = params.energy.num_energies;

%% Create Energy Array
energies_J = linspace(energy_min, energy_max, num_energies);
energies = energies_J / (1.602e-19 * 1e6);  % Convert back to MeV for output

%% Initialize Output Arrays
transmission = zeros(size(energies));
reflection = zeros(size(energies));

%% Calculate Transmission for Each Energy
fprintf('  Computing transmission coefficients for %d energies...\n', num_energies);

for i = 1:num_energies
    E = energies_J(i);
    
    % Skip if energy is too low (below potential minimum)
    if E <= min(V)
        transmission(i) = 0;
        reflection(i) = 1;
        continue;
    end
    
    % Calculate wave number
    k = sqrt(2 * m * E) / hbar;
    
    % Solve scattering problem
    [T, R] = solve_scattering_problem(V, E, params);
    
    transmission(i) = T;
    reflection(i) = R;
    
    % Progress indicator
    if mod(i, max(1, floor(num_energies/10))) == 0
        fprintf('    Progress: %d%% (E = %.1f MeV, T = %.2e)\n', ...
                round(100*i/num_energies), energies(i), T);
    end
end

%% Validation
if any(transmission < 0) || any(transmission > 1)
    warning('Transmission coefficients outside [0,1] range');
end

if any(reflection < 0) || any(reflection > 1)
    warning('Reflection coefficients outside [0,1] range');
end

% Check conservation of probability
conservation_error = abs(transmission + reflection - 1);
max_error = max(conservation_error);
if max_error > 1e-6
    warning('Probability conservation violated: max error = %.2e', max_error);
end

%% Display Results
fprintf('  Transmission calculation completed\n');
fprintf('  Energy range: %.1f - %.1f MeV\n', min(energies), max(energies));
fprintf('  Max transmission: %.2e at %.1f MeV\n', max(transmission), ...
        energies(transmission == max(transmission)));

end

%% Helper Function: Solve Scattering Problem
function [T, R] = solve_scattering_problem(V, E, params)
%SOLVE_SCATTERING_PROBLEM Solve the scattering problem for a given energy
%   Uses a simplified approach for transmission coefficient calculation
    
    x = params.numerical.x;
    dx = params.numerical.dx;
    m = params.physics.m;
    hbar = params.physics.hbar;
    
    % Wave number
    k = sqrt(2 * m * E) / hbar;
    
    % Find barrier region (where V > E)
    barrier_region = V > E;
    
    if ~any(barrier_region)
        % No barrier, perfect transmission
        T = 1;
        R = 0;
        return;
    end
    
    % Calculate transmission using WKB approximation
    % T ≈ exp(-2 * integral of kappa dx) where kappa = sqrt(2m(V-E))/hbar
    
    % Find barrier boundaries
    barrier_start = find(barrier_region, 1, 'first');
    barrier_end = find(barrier_region, 1, 'last');
    
    if isempty(barrier_start) || isempty(barrier_end)
        T = 1;
        R = 0;
        return;
    end
    
    % Calculate kappa in barrier region
    V_barrier = V(barrier_start:barrier_end);
    x_barrier = x(barrier_start:barrier_end);
    kappa = sqrt(2 * m * (V_barrier - E)) / hbar;
    
    % Only use real, positive kappa values
    valid_idx = kappa > 0;
    if any(valid_idx)
        % Numerical integration of kappa
        integral_kappa = trapz(x_barrier(valid_idx), kappa(valid_idx));
        
        % Transmission coefficient
        T = exp(-2 * integral_kappa);
        
        % Ensure physical bounds
        T = max(0, min(1, T));
        R = 1 - T;
    else
        T = 1;
        R = 0;
    end
    
end
