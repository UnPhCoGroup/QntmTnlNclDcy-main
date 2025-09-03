function main(varargin)
%MAIN Main function for quantum tunneling simulation
%
%   This is the central orchestrator that coordinates all components of the
%   quantum tunneling simulation. It initializes parameters, sets up the
%   potential, solves the Schrödinger equation, computes transmission
%   coefficients, and generates visualizations.
%
%   Inputs:
%       varargin - Optional parameter-value pairs passed to init_params
%
%   Outputs:
%       None (saves results to data/ directory and generates plots)
%
%   Workflow:
%       1. Initialize parameters using init_params
%       2. Generate potential using potential_square or potential_coulomb
%       3. Solve Schrödinger equation using solve_schrodinger
%       4. Compute transmission coefficients using compute_transmission
%       5. Visualize results using visualize_wavefunction
%       6. Save results to data/results/
%
%   Example:
%       main();  % Run with default parameters
%       main('barrier_height', 40, 'energy', 10);  % Custom parameters
%
%   Author: UniPhi Collective
%   Date: July-September 2025

%% Clear workspace and close figures
clear; clc; close all;

fprintf('Starting Quantum Tunneling Simulation...\n');
fprintf('=====================================\n');

%% Step 1: Initialize Parameters
fprintf('Step 1: Initializing parameters...\n');
params = init_params(varargin{:});

%% Step 2: Generate Potential
fprintf('Step 2: Generating potential...\n');
switch params.simulation.barrier_type
    case 'square'
        V = potential_square(params);
        fprintf('  Using square barrier potential\n');
    case 'double'
        V = potential_double_barrier(params);
        fprintf('  Using double barrier potential\n');
    case 'coulomb'
        V = potential_coulomb(params);
        fprintf('  Using Coulomb barrier potential\n');
    otherwise
        error('Unknown barrier type: %s', params.simulation.barrier_type);
end

%% Step 3: Solve Schrödinger Equation
fprintf('Step 3: Solving Schrödinger equation...\n');
tic;
[eigenvalues, eigenfunctions] = solve_schrodinger(V, params);
solve_time = toc;
fprintf('  Solver completed in %.2f seconds\n', solve_time);
fprintf('  Found %d bound states\n', length(eigenvalues));

%% Step 4: Compute Transmission Coefficients
fprintf('Step 4: Computing transmission coefficients...\n');
tic;
[transmission, reflection, energies] = compute_transmission(V, params);
transmission_time = toc;
fprintf('  Transmission calculation completed in %.2f seconds\n', transmission_time);

%% Step 5: Generate Visualizations
fprintf('Step 5: Generating visualizations...\n');
visualize_wavefunction(eigenvalues, eigenfunctions, V, params);
visualize_transmission(transmission, reflection, energies, params);

%% Step 6: Save Results
fprintf('Step 6: Saving results...\n');
save_results(eigenvalues, eigenfunctions, transmission, reflection, ...
             energies, V, params);

%% Step 7: Generate Summary Report
fprintf('Step 7: Generating summary report...\n');
generate_summary_report(eigenvalues, transmission, params);

fprintf('\nSimulation completed successfully!\n');
fprintf('Results saved to data/results/\n');
fprintf('Figures saved to figs/\n');

end

%% Helper Functions

function save_results(eigenvalues, eigenfunctions, transmission, reflection, ...
                      energies, V, params)
%SAVE_RESULTS Save simulation results to MAT files
    
    % Create results directory if it doesn't exist
    if ~exist('data/results', 'dir')
        mkdir('data/results');
    end
    
    % Save wavefunction data
    save('data/results/wavefunction_profiles.mat', ...
         'eigenvalues', 'eigenfunctions', 'V', 'params');
    
    % Save transmission data
    save('data/results/transmission_vs_energy.mat', ...
         'transmission', 'reflection', 'energies', 'params');
    
    fprintf('  Results saved to data/results/\n');
end

function visualize_transmission(transmission, reflection, energies, params)
%VISUALIZE_TRANSMISSION Create transmission coefficient plots
    
    % Create figures directory if it doesn't exist
    if ~exist('figs', 'dir')
        mkdir('figs');
    end
    
    figure('Position', [100, 100, 800, 600]);
    
    % Main transmission plot
    subplot(2,1,1);
    semilogy(energies, transmission, 'b-', 'LineWidth', 2);
    hold on;
    semilogy(energies, reflection, 'r-', 'LineWidth', 2);
    xlabel('Energy (MeV)');
    ylabel('Coefficient');
    title('Transmission and Reflection Coefficients');
    legend('Transmission', 'Reflection', 'Location', 'best');
    grid on;
    
    % Linear scale plot
    subplot(2,1,2);
    plot(energies, transmission, 'b-', 'LineWidth', 2);
    hold on;
    plot(energies, reflection, 'r-', 'LineWidth', 2);
    xlabel('Energy (MeV)');
    ylabel('Coefficient');
    title('Transmission and Reflection Coefficients (Linear Scale)');
    legend('Transmission', 'Reflection', 'Location', 'best');
    grid on;
    
    % Save figure
    saveas(gcf, 'figs/transmission_coefficients.png');
    saveas(gcf, 'figs/transmission_coefficients.fig');
    
    fprintf('  Transmission plots saved to figs/\n');
end

function generate_summary_report(eigenvalues, transmission, params)
%GENERATE_SUMMARY_REPORT Create a summary of simulation results
    
    fprintf('\n=== SIMULATION SUMMARY ===\n');
    fprintf('Barrier Type: %s\n', params.simulation.barrier_type);
    fprintf('Barrier Height: %.1f MeV\n', params.potential.barrier_height);
    fprintf('Incident Energy: %.1f MeV\n', params.energy.incident_energy);
    fprintf('Number of Bound States: %d\n', length(eigenvalues));
    
    if ~isempty(eigenvalues)
        fprintf('Ground State Energy: %.3f MeV\n', eigenvalues(1));
        if length(eigenvalues) > 1
            fprintf('First Excited State: %.3f MeV\n', eigenvalues(2));
        end
    end
    
    % Find transmission at incident energy
    [~, idx] = min(abs(transmission.energies - params.energy.incident_energy));
    T_at_incident = transmission.values(idx);
    fprintf('Transmission at %.1f MeV: %.2e\n', ...
            params.energy.incident_energy, T_at_incident);
    
    fprintf('========================\n');
end
