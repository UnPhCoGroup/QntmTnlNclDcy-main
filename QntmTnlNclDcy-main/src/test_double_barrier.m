%TEST_DOUBLE_BARRIER Test script for double barrier quantum tunneling
%
%   This script demonstrates the new double barrier functionality
%   inspired by the reference visualization.
%
%   Author: UniPhi Collective
%   Date: July-September 2025

clear; clc; close all;

fprintf('Testing Double Barrier Quantum Tunneling\n');
fprintf('=======================================\n');

%% Test 1: Basic Double Barrier Setup
fprintf('\nTest 1: Basic Double Barrier Setup\n');
params = init_params('barrier_type', 'double', 'gap_width', 4e-15, 'energy', 5);
V = potential_double_barrier(params);

% Plot the potential
figure('Position', [100, 100, 800, 400]);
x_fm = params.numerical.x * 1e15;
V_MeV = V / (1.602e-19 * 1e6);
plot(x_fm, V_MeV, 'r-', 'LineWidth', 2);
xlabel('Position (fm)');
ylabel('Potential Energy (MeV)');
title('Double Barrier Potential Structure');
grid on;
saveas(gcf, 'figs/double_barrier_potential.png');

%% Test 2: Solve Schrödinger Equation
fprintf('\nTest 2: Solving Schrödinger Equation\n');
[eigenvalues, eigenfunctions] = solve_schrodinger(V, params);
fprintf('Found %d bound states\n', length(eigenvalues));

%% Test 3: Visualize Wavefunctions
fprintf('\nTest 3: Visualizing Wavefunctions\n');
visualize_wavefunction(eigenvalues, eigenfunctions, V, params);

%% Test 4: James-style Single Wavefunction
fprintf('\nTest 4: James-style Single Wavefunction\n');
visualize_single_wavefunction(params.numerical.x, eigenfunctions, V, eigenvalues, 1, ...
                             'out_path', 'figs/double_barrier_single_wavefunction.png');

%% Test 5: Transmission Coefficients
fprintf('\nTest 5: Computing Transmission Coefficients\n');
[transmission, reflection, energies] = compute_transmission(V, params);

% Plot transmission vs energy
figure('Position', [200, 200, 800, 400]);
semilogy(energies, transmission, 'b-', 'LineWidth', 2);
hold on;
semilogy(energies, reflection, 'r-', 'LineWidth', 2);
xlabel('Energy (MeV)');
ylabel('Coefficient');
title('Double Barrier Transmission and Reflection');
legend('Transmission', 'Reflection', 'Location', 'best');
grid on;
saveas(gcf, 'figs/double_barrier_transmission.png');

%% Test 6: Quick Animation (Short Version)
fprintf('\nTest 6: Creating Quick Animation\n');
animate_double_barrier_tunneling(params, 'energy', 5, 'duration', 2, 'fps', 15);

%% Summary
fprintf('\n=== TEST SUMMARY ===\n');
fprintf('✅ Double barrier potential created\n');
fprintf('✅ Schrödinger equation solved\n');
fprintf('✅ Wavefunctions visualized\n');
fprintf('✅ Transmission coefficients computed\n');
fprintf('✅ Animation created\n');
fprintf('\nFiles created:\n');
fprintf('  - figs/double_barrier_potential.png\n');
fprintf('  - figs/double_barrier_single_wavefunction.png\n');
fprintf('  - figs/double_barrier_transmission.png\n');
fprintf('  - figs/double_barrier_tunneling.mp4\n');
fprintf('  - Various wavefunction analysis plots\n');
fprintf('\nDouble barrier testing completed successfully!\n');
