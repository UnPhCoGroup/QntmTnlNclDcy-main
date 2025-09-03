function visualize_single_wavefunction(x, psi, V, E, index, varargin)
%VISUALIZE_SINGLE_WAVEFUNCTION Create James-style single wavefunction visualization
%
%   This function creates a clean, professional visualization of a single
%   wavefunction with scaled potential, inspired by James's Python implementation.
%   It scales the potential to probability density units for better visualization.
%
%   Inputs:
%       x - Position array [m]
%       psi - Wavefunction matrix [num_points, num_states] or single wavefunction
%       V - Potential energy array [J]
%       E - Energy eigenvalues [J] or single energy value
%       index - Index of wavefunction to plot (if psi is matrix)
%       varargin - Optional parameter-value pairs:
%                  'out_path' - Output file path (default: 'figs/single_wavefunction.png')
%                  'title' - Plot title (default: auto-generated)
%                  'units' - Display units: 'si' or 'nuclear' (default: 'nuclear')
%
%   Outputs:
%       None (creates and saves figure)
%
%   Example:
%       params = init_params();
%       V = potential_square(params);
%       [E, psi] = solve_schrodinger(V, params);
%       visualize_single_wavefunction(params.numerical.x, psi, V, E, 0);
%
%   Author: James (UniPhi Collective) - MATLAB adaptation
%   Date: July-September 2025

%% Parse Input Arguments
p = inputParser;
addParameter(p, 'out_path', 'figs/single_wavefunction.png', @ischar);
addParameter(p, 'title', '', @ischar);
addParameter(p, 'units', 'nuclear', @(x) ismember(x, {'si', 'nuclear'}));
parse(p, varargin{:});

out_path = p.Results.out_path;
title_str = p.Results.title;
units = p.Results.units;

%% Validate and Process Inputs
x = x(:);  % Ensure column vector
V = V(:);  % Ensure column vector

% Handle wavefunction input
if size(psi, 2) > 1
    % Matrix input - select specific state
    if index < 1 || index > size(psi, 2)
        error('Index %d out of range [1, %d]', index, size(psi, 2));
    end
    psi_vec = psi(:, index);
    E_value = E(index);
else
    % Single wavefunction
    psi_vec = psi(:);
    E_value = E;
end

%% Normalize Wavefunction
norm = trapz(x, abs(psi_vec).^2);
if norm <= 0
    norm = 1e-12;
end
psi_norm = psi_vec / sqrt(norm);
prob_density = abs(psi_norm).^2;

%% Convert Units
if strcmp(units, 'nuclear')
    % Convert to nuclear units (fm, MeV)
    x_display = x * 1e15;  % Convert to fm
    V_display = V / (1.602e-19 * 1e6);  % Convert to MeV
    E_display = E_value / (1.602e-19 * 1e6);  % Convert to MeV
    x_label = 'Position (fm)';
    y_label = 'Probability Density (arb. units)';
else
    % Keep SI units
    x_display = x;
    V_display = V;
    E_display = E_value;
    x_label = 'Position (m)';
    y_label = 'Probability Density (m^{-1})';
end

%% Scale Potential to Probability Units
vmax = max(abs(V_display));
pmax = max(prob_density);
if vmax > 0 && pmax > 0
    V_scaled = V_display * (pmax / vmax);
    E_scaled = E_display * (pmax / vmax);
else
    V_scaled = zeros(size(V_display));
    E_scaled = 0;
end

%% Create Figure
fig = figure('Position', [100, 100, 800, 500]);

% Plot probability density
plot(x_display, prob_density, 'b-', 'LineWidth', 2.5, ...
     'DisplayName', '|ψ(x)|²');
hold on;

% Plot scaled potential
plot(x_display, V_scaled, 'k--', 'LineWidth', 1.8, ...
     'DisplayName', 'Scaled V(x)');

% Plot scaled energy level
yline(E_scaled, 'r--', 'LineWidth', 1.5, ...
      'DisplayName', 'Scaled E');

%% Formatting
xlabel(x_label, 'FontSize', 12);
ylabel(y_label, 'FontSize', 12);

% Auto-generate title if not provided
if isempty(title_str)
    if strcmp(units, 'nuclear')
        title_str = sprintf('Wavefunction State %d (E = %.2f MeV)', ...
                           index-1, E_display);
    else
        title_str = sprintf('Wavefunction State %d (E = %.2e J)', ...
                           index-1, E_display);
    end
end
title(title_str, 'FontSize', 14, 'FontWeight', 'bold');

% Grid and legend
grid on;
grid minor;
legend('Location', 'best', 'FontSize', 11);

% Improve appearance
set(gca, 'FontSize', 11);
set(gca, 'LineWidth', 1.2);

%% Save Figure
% Create output directory if it doesn't exist
[out_dir, ~, ~] = fileparts(out_path);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Save in multiple formats
[~, name, ~] = fileparts(out_path);
saveas(fig, out_path);
saveas(fig, strrep(out_path, '.png', '.fig'));

% Close figure to save memory
close(fig);

%% Display Results
fprintf('Single wavefunction visualization saved:\n');
fprintf('  - %s\n', out_path);
fprintf('  - %s\n', strrep(out_path, '.png', '.fig'));
if strcmp(units, 'nuclear')
    fprintf('  - State %d, Energy = %.2f MeV\n', index-1, E_display);
else
    fprintf('  - State %d, Energy = %.2e J\n', index-1, E_display);
end

end
