function [eigenvalues, eigenfunctions] = solve_schrodinger(V, params)
%SOLVE_SCHRODINGER Solve the 1D time-independent Schrödinger equation
%
%   This function solves the time-independent Schrödinger equation using
%   finite difference methods to find eigenvalues and eigenfunctions for
%   a given potential V(x).
%
%   Inputs:
%       V - Potential energy array [J] at each grid point
%       params - Structure containing simulation parameters
%                Required fields: numerical.x, numerical.dx, physics.m, physics.hbar
%
%   Outputs:
%       eigenvalues - Array of energy eigenvalues [J] (sorted in ascending order)
%       eigenfunctions - Matrix where each column is an eigenfunction
%                        Size: [num_points, num_eigenvalues]
%
%   Method:
%       Uses finite difference discretization of the Hamiltonian:
%       H = -ℏ²/(2m) * d²/dx² + V(x)
%       
%       The second derivative is approximated using central differences:
%       d²ψ/dx² ≈ (ψ[i+1] - 2ψ[i] + ψ[i-1]) / dx²
%
%   Example:
%       params = init_params();
%       V = potential_square(params);
%       [E, psi] = solve_schrodinger(V, params);
%
%   Author: Alex (UniPhi Collective)
%   Date: July-September 2025

%% Validate Input
if ~isfield(params, 'numerical') || ~isfield(params.numerical, 'x')
    error('params must contain numerical.x field');
end

if ~isfield(params, 'physics') || ~isfield(params.physics, 'm')
    error('params must contain physics.m field');
end

if length(V) ~= length(params.numerical.x)
    error('Potential V must have same length as params.numerical.x');
end

%% Extract Parameters
x = params.numerical.x;
dx = params.numerical.dx;
m = params.physics.m;
hbar = params.physics.hbar;
num_points = length(x);

%% Construct Hamiltonian Matrix
fprintf('  Constructing Hamiltonian matrix...\n');

% Kinetic energy coefficient
T_coeff = -hbar^2 / (2 * m * dx^2);

% Initialize Hamiltonian matrix (sparse for efficiency)
H = zeros(num_points, num_points);

% Fill kinetic energy terms (tridiagonal structure)
for i = 2:num_points-1
    H(i, i-1) = T_coeff;      % Lower diagonal
    H(i, i) = -2 * T_coeff;   % Main diagonal
    H(i, i+1) = T_coeff;      % Upper diagonal
end

% Add potential energy terms to diagonal
for i = 1:num_points
    H(i, i) = H(i, i) + V(i);
end

%% Apply Boundary Conditions
% For bound states, we typically use zero boundary conditions
% ψ(x_min) = ψ(x_max) = 0
% This is already satisfied by our finite difference scheme

% Alternative: For scattering states, we might use different boundary conditions
% For now, we'll use zero boundary conditions for bound states

%% Solve Eigenvalue Problem
fprintf('  Solving eigenvalue problem...\n');
tic;

% Use MATLAB's built-in eigensolver
% We only need the lowest few eigenvalues and eigenvectors
num_eigenstates = min(10, num_points);  % Limit to first 10 states

[eigenfunctions, eigenvals] = eigs(H, num_eigenstates, 'smallestabs');
eigenvalues = diag(eigenvals);

eigenvalue_time = toc;
fprintf('  Eigenvalue solution completed in %.2f seconds\n', eigenvalue_time);

%% Sort Results
[eigenvalues, sort_idx] = sort(eigenvalues);
eigenfunctions = eigenfunctions(:, sort_idx);

%% Filter Bound States
% Only keep states with negative energy (bound states)
bound_states = eigenvalues < 0;
eigenvalues = eigenvalues(bound_states);
eigenfunctions = eigenfunctions(:, bound_states);

%% Normalize Eigenfunctions
fprintf('  Normalizing eigenfunctions...\n');
for i = 1:size(eigenfunctions, 2)
    % Normalize using trapezoidal rule
    norm_factor = sqrt(trapz(x, abs(eigenfunctions(:, i)).^2));
    if norm_factor > 0
        eigenfunctions(:, i) = eigenfunctions(:, i) / norm_factor;
    end
end

%% Validation
if any(isnan(eigenvalues)) || any(isinf(eigenvalues))
    error('Eigenvalues contain NaN or Inf values');
end

if any(any(isnan(eigenfunctions))) || any(any(isinf(eigenfunctions)))
    error('Eigenfunctions contain NaN or Inf values');
end

%% Display Results
fprintf('  Found %d bound states\n', length(eigenvalues));
if ~isempty(eigenvalues)
    fprintf('  Ground state energy: %.3f MeV\n', ...
            eigenvalues(1) / (1.602e-19 * 1e6));
    if length(eigenvalues) > 1
        fprintf('  First excited state: %.3f MeV\n', ...
                eigenvalues(2) / (1.602e-19 * 1e6));
    end
end

%% Optional: Check Orthogonality
if size(eigenfunctions, 2) > 1
    fprintf('  Checking orthogonality...\n');
    overlap_matrix = eigenfunctions' * eigenfunctions;
    max_off_diagonal = max(max(abs(overlap_matrix - eye(size(overlap_matrix)))));
    fprintf('  Maximum off-diagonal overlap: %.2e\n', max_off_diagonal);
    
    if max_off_diagonal > 1e-6
        warning('Eigenfunctions may not be properly orthogonal');
    end
end

end
