function visualize_wavefunction(eigenvalues, eigenfunctions, V, params)
%VISUALIZE_WAVEFUNCTION Create comprehensive wavefunction visualizations
%
%   This function generates multiple plots to visualize the quantum
%   mechanical wavefunctions, probability densities, and potential energy
%   landscape. It creates both static plots and saves them to the figs/
%   directory.
%
%   Inputs:
%       eigenvalues - Array of energy eigenvalues [J]
%       eigenfunctions - Matrix of eigenfunctions [num_points, num_states]
%       V - Potential energy array [J] at each grid point
%       params - Structure containing simulation parameters
%
%   Outputs:
%       None (creates and saves figures)
%
%   Generated Plots:
%       1. Potential energy and energy levels
%       2. Wavefunctions and probability densities
%       3. Combined potential and wavefunction plot
%       4. Energy level diagram
%
%   Example:
%       params = init_params();
%       V = potential_square(params);
%       [E, psi] = solve_schrodinger(V, params);
%       visualize_wavefunction(E, psi, V, params);
%
%   Author: James (UniPhi Collective)
%   Date: July-September 2025

%% Validate Input
if nargin < 4
    error('All four arguments required: eigenvalues, eigenfunctions, V, params');
end

if isempty(eigenvalues) || isempty(eigenfunctions)
    warning('No eigenstates to visualize');
    return;
end

%% Extract Parameters
x = params.numerical.x;
x_fm = x * 1e15;  % Convert to femtometers for display
V_MeV = V / (1.602e-19 * 1e6);  % Convert to MeV for display
E_MeV = eigenvalues / (1.602e-19 * 1e6);  % Convert to MeV for display

%% Helper Function: Scale Potential to Probability Units
function [V_scaled, E_scaled] = scale_to_prob_units(V, prob, E_value)
    vmax = max(abs(V));
    pmax = max(prob);
    if vmax > 0 && pmax > 0
        V_scaled = V * (pmax / vmax);
        E_scaled = E_value * (pmax / vmax);
    else
        V_scaled = zeros(size(V));
        E_scaled = 0;
    end
end

%% Create Figures Directory
if ~exist('figs', 'dir')
    mkdir('figs');
end

%% Plot 1: Potential Energy and Energy Levels
figure('Position', [100, 100, 1000, 600]);
subplot(2,2,1);

% Plot potential
plot(x_fm, V_MeV, 'k-', 'LineWidth', 2);
hold on;

% Plot energy levels
for i = 1:length(E_MeV)
    yline(E_MeV(i), 'r--', 'LineWidth', 1.5);
    text(max(x_fm) * 0.8, E_MeV(i), sprintf('E_%d = %.2f MeV', i-1, E_MeV(i)), ...
         'FontSize', 8, 'BackgroundColor', 'white');
end

xlabel('Position (fm)');
ylabel('Energy (MeV)');
title('Potential Energy and Bound State Energies');
grid on;
legend('Potential V(x)', 'Energy Levels', 'Location', 'best');

%% Plot 2: Wavefunctions
subplot(2,2,2);

% Plot first few wavefunctions
colors = {'b-', 'r-', 'g-', 'm-', 'c-'};
for i = 1:min(5, size(eigenfunctions, 2))
    psi_shifted = eigenfunctions(:, i) + E_MeV(i);  % Shift for visibility
    plot(x_fm, psi_shifted, colors{mod(i-1, length(colors)) + 1}, ...
         'LineWidth', 1.5, 'DisplayName', sprintf('ψ_%d', i-1));
    hold on;
end

xlabel('Position (fm)');
ylabel('Energy (MeV)');
title('Wavefunctions (Shifted for Visibility)');
grid on;
legend('Location', 'best');

%% Plot 3: Probability Densities
subplot(2,2,3);

% Plot probability densities
for i = 1:min(5, size(eigenfunctions, 2))
    prob_density = abs(eigenfunctions(:, i)).^2;
    prob_density = prob_density / max(prob_density) * 0.5;  % Normalize for visibility
    prob_density = prob_density + E_MeV(i);  % Shift for visibility
    
    plot(x_fm, prob_density, colors{mod(i-1, length(colors)) + 1}, ...
         'LineWidth', 1.5, 'DisplayName', sprintf('|ψ_%d|²', i-1));
    hold on;
end

xlabel('Position (fm)');
ylabel('Energy (MeV)');
title('Probability Densities (Normalized and Shifted)');
grid on;
legend('Location', 'best');

%% Plot 4: Combined View
subplot(2,2,4);

% Plot potential
plot(x_fm, V_MeV, 'k-', 'LineWidth', 2);
hold on;

% Plot energy levels
for i = 1:length(E_MeV)
    yline(E_MeV(i), 'r--', 'LineWidth', 1);
end

% Plot ground state wavefunction
if size(eigenfunctions, 2) >= 1
    psi_0 = eigenfunctions(:, 1);
    psi_0_shifted = psi_0 / max(abs(psi_0)) * 2 + E_MeV(1);  % Scale and shift
    plot(x_fm, psi_0_shifted, 'b-', 'LineWidth', 1.5);
end

xlabel('Position (fm)');
ylabel('Energy (MeV)');
title('Potential, Energy Levels, and Ground State');
grid on;
legend('Potential', 'Energy Levels', 'Ground State ψ₀', 'Location', 'best');

% Save combined figure
sgtitle(sprintf('Quantum Tunneling Simulation: %s Barrier', ...
                params.simulation.barrier_type), 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'figs/wavefunction_analysis.png');
saveas(gcf, 'figs/wavefunction_analysis.fig');

%% Create Detailed Wavefunction Plot
figure('Position', [200, 200, 1200, 800]);

% Plot potential
yyaxis left;
plot(x_fm, V_MeV, 'k-', 'LineWidth', 3);
ylabel('Potential Energy (MeV)');
ylim([min(V_MeV) - 5, max(V_MeV) + 5]);

% Plot wavefunctions
yyaxis right;
colors = {'b-', 'r-', 'g-', 'm-', 'c-', 'y-'};
for i = 1:min(6, size(eigenfunctions, 2))
    psi_shifted = eigenfunctions(:, i) * 3 + E_MeV(i);  % Scale and shift
    plot(x_fm, psi_shifted, colors{i}, 'LineWidth', 2, ...
         'DisplayName', sprintf('ψ_%d (E = %.2f MeV)', i-1, E_MeV(i)));
    hold on;
end

xlabel('Position (fm)');
ylabel('Energy (MeV)');
title(sprintf('Detailed Wavefunction Analysis - %s Barrier', ...
              params.simulation.barrier_type), 'FontSize', 16, 'FontWeight', 'bold');
grid on;
legend('Location', 'best', 'FontSize', 10);

% Add energy level markers
for i = 1:length(E_MeV)
    text(max(x_fm) * 0.02, E_MeV(i), sprintf('E_%d', i-1), ...
         'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'white');
end

saveas(gcf, 'figs/detailed_wavefunctions.png');
saveas(gcf, 'figs/detailed_wavefunctions.fig');

%% Create Energy Level Diagram
figure('Position', [300, 300, 600, 800]);

% Create energy level diagram
y_positions = E_MeV;
x_positions = zeros(size(y_positions));

% Plot energy levels as horizontal lines
for i = 1:length(E_MeV)
    plot([-1, 1], [y_positions(i), y_positions(i)], 'r-', 'LineWidth', 3);
    hold on;
    
    % Add labels
    text(0.1, y_positions(i), sprintf('n = %d', i-1), ...
         'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    text(-0.1, y_positions(i), sprintf('%.3f MeV', E_MeV(i)), ...
         'FontSize', 10, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end

% Add potential well representation
well_bottom = min(V_MeV);
well_top = max(V_MeV);
plot([-0.5, 0.5], [well_bottom, well_bottom], 'k-', 'LineWidth', 4);
plot([-0.5, 0.5], [well_top, well_top], 'k-', 'LineWidth', 4);
plot([-0.5, -0.5], [well_bottom, well_top], 'k-', 'LineWidth', 4);
plot([0.5, 0.5], [well_bottom, well_top], 'k-', 'LineWidth', 4);

xlim([-1.5, 1.5]);
ylim([well_bottom - 2, well_top + 2]);
xlabel('Quantum Number');
ylabel('Energy (MeV)');
title('Energy Level Diagram', 'FontSize', 16, 'FontWeight', 'bold');
grid on;

saveas(gcf, 'figs/energy_levels.png');
saveas(gcf, 'figs/energy_levels.fig');

%% Create James-style Scaled Visualization
figure('Position', [400, 400, 1000, 600]);

% Plot first few states with scaled potential
colors = {'b-', 'r-', 'g-', 'm-', 'c-'};
for i = 1:min(5, size(eigenfunctions, 2))
    subplot(2, 3, i);
    
    % Get probability density
    prob_density = abs(eigenfunctions(:, i)).^2;
    
    % Scale potential to probability units
    [V_scaled, E_scaled] = scale_to_prob_units(V_MeV, prob_density, E_MeV(i));
    
    % Plot
    plot(x_fm, prob_density, colors{mod(i-1, length(colors)) + 1}, 'LineWidth', 2);
    hold on;
    plot(x_fm, V_scaled, 'k--', 'LineWidth', 1.5);
    yline(E_scaled, 'r--', 'LineWidth', 1.2);
    
    xlabel('Position (fm)');
    ylabel('Probability Density');
    title(sprintf('State %d (E = %.2f MeV)', i-1, E_MeV(i)));
    legend('|ψ|²', 'Scaled V(x)', 'Scaled E', 'Location', 'best');
    grid on;
end

% Summary subplot
subplot(2, 3, 6);
% Plot all probability densities together
for i = 1:min(5, size(eigenfunctions, 2))
    prob_density = abs(eigenfunctions(:, i)).^2;
    plot(x_fm, prob_density, colors{mod(i-1, length(colors)) + 1}, ...
         'LineWidth', 1.5, 'DisplayName', sprintf('State %d', i-1));
    hold on;
end
xlabel('Position (fm)');
ylabel('Probability Density');
title('All Bound States');
legend('Location', 'best');
grid on;

sgtitle('James-Style Scaled Wavefunction Visualization', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'figs/james_style_wavefunctions.png');
saveas(gcf, 'figs/james_style_wavefunctions.fig');

%% Display Summary
fprintf('  Wavefunction visualizations created:\n');
fprintf('    - Combined analysis plot: figs/wavefunction_analysis.png\n');
fprintf('    - Detailed wavefunctions: figs/detailed_wavefunctions.png\n');
fprintf('    - Energy level diagram: figs/energy_levels.png\n');
fprintf('    - James-style scaled plots: figs/james_style_wavefunctions.png\n');
fprintf('    - Number of bound states visualized: %d\n', size(eigenfunctions, 2));

end
