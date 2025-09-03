function results = sweep_parameters(param_name, param_values, varargin)
%SWEEP_PARAMETERS Perform parameter sweep analysis for quantum tunneling
%
%   This function systematically varies a specified parameter and analyzes
%   how it affects the quantum tunneling simulation results. It's useful
%   for understanding the sensitivity of the system to different parameters
%   and for finding optimal parameter values.
%
%   Inputs:
%       param_name - Name of parameter to sweep (string)
%                    Supported: 'barrier_height', 'barrier_width', 'energy',
%                              'well_width', 'particle_mass'
%       param_values - Array of parameter values to test
%       varargin - Optional parameter-value pairs for fixed parameters
%
%   Outputs:
%       results - Structure containing sweep results:
%                 .param_name - Name of swept parameter
%                 .param_values - Array of parameter values tested
%                 .eigenvalues - Matrix of eigenvalues [num_values, num_states]
%                 .transmission_max - Maximum transmission coefficient for each value
%                 .transmission_at_incident - Transmission at incident energy
%                 .ground_state_energy - Ground state energy for each value
%                 .num_bound_states - Number of bound states for each value
%
%   Example:
%       % Sweep barrier height from 10 to 50 MeV
%       barrier_heights = 10:5:50;
%       results = sweep_parameters('barrier_height', barrier_heights);
%       
%       % Sweep energy with fixed barrier height
%       energies = 1:0.5:20;
%       results = sweep_parameters('energy', energies, 'barrier_height', 30);
%
%   Author: Alex (UniPhi Collective)
%   Date: July-September 2025

%% Validate Input
if nargin < 2
    error('At least two arguments required: param_name and param_values');
end

if ~ischar(param_name) && ~isstring(param_name)
    error('param_name must be a string');
end

if ~isnumeric(param_values) || length(param_values) < 2
    error('param_values must be a numeric array with at least 2 elements');
end

%% Initialize Results Structure
results.param_name = param_name;
results.param_values = param_values;
num_values = length(param_values);

% Pre-allocate result arrays
results.eigenvalues = [];
results.transmission_max = zeros(num_values, 1);
results.transmission_at_incident = zeros(num_values, 1);
results.ground_state_energy = zeros(num_values, 1);
results.num_bound_states = zeros(num_values, 1);
results.computation_time = zeros(num_values, 1);

%% Display Progress
fprintf('Starting parameter sweep: %s\n', param_name);
fprintf('Number of parameter values: %d\n', num_values);
fprintf('Parameter range: %.3f to %.3f\n', min(param_values), max(param_values));
fprintf('=====================================\n');

%% Perform Parameter Sweep
for i = 1:num_values
    fprintf('Step %d/%d: %s = %.3f\n', i, num_values, param_name, param_values(i));
    
    tic;
    
    try
        %% Set Up Parameters
        % Create parameter structure with current sweep value
        sweep_params = varargin;
        sweep_params{end+1} = param_name;
        sweep_params{end+1} = param_values(i);
        
        % Initialize parameters
        params = init_params(sweep_params{:});
        
        %% Generate Potential
        switch params.simulation.barrier_type
            case 'square'
                V = potential_square(params);
            case 'double'
                V = potential_double_barrier(params);
            case 'coulomb'
                V = potential_coulomb(params);
            otherwise
                error('Unknown barrier type: %s', params.simulation.barrier_type);
        end
        
        %% Solve Schrödinger Equation
        [eigenvalues, eigenfunctions] = solve_schrodinger(V, params);
        
        %% Compute Transmission Coefficients
        [transmission, reflection, energies] = compute_transmission(V, params);
        
        %% Store Results
        if i == 1
            % Initialize eigenvalue matrix on first iteration
            max_states = 10;  % Maximum number of states to track
            results.eigenvalues = zeros(num_values, max_states);
        end
        
        % Store eigenvalues (pad with NaN if fewer states)
        num_states = min(length(eigenvalues), max_states);
        results.eigenvalues(i, 1:num_states) = eigenvalues(1:num_states);
        if num_states < max_states
            results.eigenvalues(i, num_states+1:end) = NaN;
        end
        
        % Store other results
        results.ground_state_energy(i) = eigenvalues(1) / (1.602e-19 * 1e6);  % Convert to MeV
        results.num_bound_states(i) = length(eigenvalues);
        results.transmission_max(i) = max(transmission);
        
        % Find transmission at incident energy
        [~, idx] = min(abs(energies - params.energy.incident_energy));
        results.transmission_at_incident(i) = transmission(idx);
        
        results.computation_time(i) = toc;
        
        fprintf('  Completed in %.2f seconds\n', results.computation_time(i));
        fprintf('  Ground state: %.3f MeV, Bound states: %d\n', ...
                results.ground_state_energy(i), results.num_bound_states(i));
        fprintf('  Max transmission: %.2e\n', results.transmission_max(i));
        
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        results.computation_time(i) = toc;
        
        % Fill with NaN for failed calculations
        if i == 1
            max_states = 10;
            results.eigenvalues = zeros(num_values, max_states);
        end
        results.eigenvalues(i, :) = NaN;
        results.ground_state_energy(i) = NaN;
        results.num_bound_states(i) = NaN;
        results.transmission_max(i) = NaN;
        results.transmission_at_incident(i) = NaN;
    end
    
    fprintf('\n');
end

%% Generate Summary Statistics
fprintf('Parameter sweep completed!\n');
fprintf('=====================================\n');
fprintf('Total computation time: %.2f seconds\n', sum(results.computation_time));
fprintf('Average time per calculation: %.2f seconds\n', mean(results.computation_time));
fprintf('Successful calculations: %d/%d\n', sum(~isnan(results.ground_state_energy)), num_values);

%% Create Visualization
visualize_parameter_sweep(results);

%% Save Results
if ~exist('data/results', 'dir')
    mkdir('data/results');
end

save('data/results/parameter_sweep.mat', 'results');
fprintf('Results saved to data/results/parameter_sweep.mat\n');

end

%% Helper Function: Visualize Parameter Sweep Results
function visualize_parameter_sweep(results)
%VISUALIZE_PARAMETER_SWEEP Create plots for parameter sweep analysis
    
    if ~exist('figs', 'dir')
        mkdir('figs');
    end
    
    param_name = results.param_name;
    param_values = results.param_values;
    
    %% Create Multi-Panel Figure
    figure('Position', [100, 100, 1200, 800]);
    
    % Plot 1: Ground State Energy
    subplot(2,3,1);
    plot(param_values, results.ground_state_energy, 'bo-', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel(param_name);
    ylabel('Ground State Energy (MeV)');
    title('Ground State Energy vs Parameter');
    grid on;
    
    % Plot 2: Number of Bound States
    subplot(2,3,2);
    plot(param_values, results.num_bound_states, 'ro-', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel(param_name);
    ylabel('Number of Bound States');
    title('Bound States vs Parameter');
    grid on;
    
    % Plot 3: Maximum Transmission
    subplot(2,3,3);
    semilogy(param_values, results.transmission_max, 'go-', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel(param_name);
    ylabel('Maximum Transmission Coefficient');
    title('Max Transmission vs Parameter');
    grid on;
    
    % Plot 4: Transmission at Incident Energy
    subplot(2,3,4);
    semilogy(param_values, results.transmission_at_incident, 'mo-', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel(param_name);
    ylabel('Transmission at Incident Energy');
    title('Transmission at Incident Energy vs Parameter');
    grid on;
    
    % Plot 5: Computation Time
    subplot(2,3,5);
    plot(param_values, results.computation_time, 'co-', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel(param_name);
    ylabel('Computation Time (seconds)');
    title('Computation Time vs Parameter');
    grid on;
    
    % Plot 6: Energy Level Diagram
    subplot(2,3,6);
    hold on;
    colors = {'b', 'r', 'g', 'm', 'c', 'y', 'k'};
    for i = 1:min(7, size(results.eigenvalues, 2))
        valid_idx = ~isnan(results.eigenvalues(:, i));
        if any(valid_idx)
            plot(param_values(valid_idx), results.eigenvalues(valid_idx, i) / (1.602e-19 * 1e6), ...
                 [colors{i} 'o-'], 'LineWidth', 1.5, 'MarkerSize', 4, ...
                 'DisplayName', sprintf('State %d', i-1));
        end
    end
    xlabel(param_name);
    ylabel('Energy (MeV)');
    title('Energy Levels vs Parameter');
    legend('Location', 'best');
    grid on;
    
    % Overall title
    sgtitle(sprintf('Parameter Sweep Analysis: %s', param_name), ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    % Save figure
    saveas(gcf, 'figs/parameter_sweep_analysis.png');
    saveas(gcf, 'figs/parameter_sweep_analysis.fig');
    
    fprintf('Parameter sweep visualization saved to figs/parameter_sweep_analysis.png\n');
    
end
