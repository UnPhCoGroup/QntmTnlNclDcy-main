function animate_tunneling(V, params, varargin)
%ANIMATE_TUNNELING Create animated visualization of quantum tunneling
%
%   This function creates an animated visualization showing the time evolution
%   of a wave packet as it encounters and tunnels through a potential barrier.
%   The animation demonstrates the quantum mechanical nature of tunneling.
%
%   Inputs:
%       V - Potential energy array [J] at each grid point
%       params - Structure containing simulation parameters
%       varargin - Optional parameter-value pairs:
%                  'energy' - Energy of the wave packet [MeV] (default: 5)
%                  'width' - Width of the wave packet (default: 2e-15)
%                  'position' - Initial position of wave packet [m] (default: -10e-15)
%                  'duration' - Animation duration in seconds (default: 5)
%                  'fps' - Frames per second (default: 30)
%
%   Outputs:
%       None (creates and saves animation)
%
%   Method:
%       Creates a Gaussian wave packet and propagates it through the potential
%       using the time-dependent Schrödinger equation. The animation shows
%       both the real and imaginary parts of the wavefunction, as well as
%       the probability density.
%
%   Example:
%       params = init_params();
%       V = potential_square(params);
%       animate_tunneling(V, params, 'energy', 10, 'duration', 3);
%
%   Author: James (UniPhi Collective)
%   Date: July-September 2025

%% Parse Input Arguments
p = inputParser;
addParameter(p, 'energy', 5, @(x) x > 0);  % MeV
addParameter(p, 'width', 2e-15, @(x) x > 0);  % m
addParameter(p, 'position', -10e-15, @isnumeric);  % m
addParameter(p, 'duration', 5, @(x) x > 0);  % seconds
addParameter(p, 'fps', 30, @(x) x > 0);
parse(p, varargin{:});

energy_MeV = p.Results.energy;
packet_width = p.Results.width;
initial_position = p.Results.position;
duration = p.Results.duration;
fps = p.Results.fps;

%% Extract Parameters
x = params.numerical.x;
dx = params.numerical.dx;
m = params.physics.m;
hbar = params.physics.hbar;
num_points = length(x);

% Convert energy to Joules
E = energy_MeV * 1.602e-19 * 1e6;

%% Create Initial Wave Packet
fprintf('  Creating initial wave packet...\n');

% Gaussian wave packet
sigma = packet_width;  % Width parameter
x0 = initial_position;  % Initial position

% Calculate wave number
k0 = sqrt(2 * m * E) / hbar;

% Create initial wavefunction
psi_initial = exp(-((x - x0).^2) / (2 * sigma^2)) .* exp(1i * k0 * x);

% Normalize
psi_initial = psi_initial / sqrt(trapz(x, abs(psi_initial).^2));

%% Time Evolution Parameters
dt = 1e-21;  % Time step [s]
num_frames = round(duration * fps);
time_points = 0:dt:(num_frames-1)*dt;

%% Construct Time Evolution Operator
fprintf('  Constructing time evolution operator...\n');

% Kinetic energy operator
T_coeff = -hbar^2 / (2 * m * dx^2);

% Create Hamiltonian matrix
H = zeros(num_points, num_points);
for i = 2:num_points-1
    H(i, i-1) = T_coeff;
    H(i, i) = -2 * T_coeff + V(i);
    H(i, i+1) = T_coeff;
end

% Time evolution operator: U = exp(-i*H*dt/ℏ)
U = expm(-1i * H * dt / hbar);

%% Create Animation
fprintf('  Creating animation (%d frames)...\n', num_frames);

% Create figures directory
if ~exist('figs', 'dir')
    mkdir('figs');
end

% Initialize wavefunction
psi = psi_initial;

% Create figure
fig = figure('Position', [100, 100, 1200, 800]);

% Pre-allocate arrays for efficiency
prob_density = zeros(num_points, num_frames);
real_part = zeros(num_points, num_frames);
imag_part = zeros(num_points, num_frames);

% Time evolution loop
for frame = 1:num_frames
    % Store current state
    prob_density(:, frame) = abs(psi).^2;
    real_part(:, frame) = real(psi);
    imag_part(:, frame) = imag(psi);
    
    % Evolve wavefunction
    psi = U * psi;
    
    % Progress indicator
    if mod(frame, max(1, floor(num_frames/10))) == 0
        fprintf('    Progress: %d%%\n', round(100*frame/num_frames));
    end
end

%% Create Animation Plots
fprintf('  Generating animation frames...\n');

% Convert to display units
x_fm = x * 1e15;
V_MeV = V / (1.602e-19 * 1e6);

% Create video writer
video_file = 'figs/tunneling_animation.mp4';
v = VideoWriter(video_file, 'MPEG-4');
v.FrameRate = fps;
open(v);

for frame = 1:num_frames
    clf;
    
    % Subplot 1: Probability density
    subplot(2,2,1);
    plot(x_fm, V_MeV, 'k-', 'LineWidth', 2);
    hold on;
    plot(x_fm, prob_density(:, frame) * max(V_MeV), 'b-', 'LineWidth', 2);
    xlabel('Position (fm)');
    ylabel('Energy (MeV)');
    title('Probability Density |ψ|²');
    legend('Potential', 'Probability Density', 'Location', 'best');
    grid on;
    ylim([min(V_MeV), max(V_MeV)]);
    
    % Subplot 2: Real part
    subplot(2,2,2);
    plot(x_fm, V_MeV, 'k-', 'LineWidth', 2);
    hold on;
    plot(x_fm, real_part(:, frame) * max(V_MeV) + max(V_MeV)/2, 'r-', 'LineWidth', 2);
    xlabel('Position (fm)');
    ylabel('Energy (MeV)');
    title('Real Part Re(ψ)');
    legend('Potential', 'Real Part', 'Location', 'best');
    grid on;
    ylim([min(V_MeV), max(V_MeV)]);
    
    % Subplot 3: Imaginary part
    subplot(2,2,3);
    plot(x_fm, V_MeV, 'k-', 'LineWidth', 2);
    hold on;
    plot(x_fm, imag_part(:, frame) * max(V_MeV) + max(V_MeV)/2, 'g-', 'LineWidth', 2);
    xlabel('Position (fm)');
    ylabel('Energy (MeV)');
    title('Imaginary Part Im(ψ)');
    legend('Potential', 'Imaginary Part', 'Location', 'best');
    grid on;
    ylim([min(V_MeV), max(V_MeV)]);
    
    % Subplot 4: Combined view
    subplot(2,2,4);
    plot(x_fm, V_MeV, 'k-', 'LineWidth', 2);
    hold on;
    plot(x_fm, prob_density(:, frame) * max(V_MeV), 'b-', 'LineWidth', 2);
    plot(x_fm, real_part(:, frame) * max(V_MeV)/2 + max(V_MeV)/2, 'r--', 'LineWidth', 1);
    plot(x_fm, imag_part(:, frame) * max(V_MeV)/2 + max(V_MeV)/2, 'g--', 'LineWidth', 1);
    xlabel('Position (fm)');
    ylabel('Energy (MeV)');
    title('Combined View');
    legend('Potential', '|ψ|²', 'Re(ψ)', 'Im(ψ)', 'Location', 'best');
    grid on;
    ylim([min(V_MeV), max(V_MeV)]);
    
    % Add time information
    current_time = time_points(frame) * 1e21;  % Convert to zeptoseconds
    sgtitle(sprintf('Quantum Tunneling Animation - t = %.1f zs (E = %.1f MeV)', ...
                    current_time, energy_MeV), 'FontSize', 14, 'FontWeight', 'bold');
    
    % Write frame to video
    writeVideo(v, getframe(fig));
end

close(v);

%% Create Static Summary Plot
figure('Position', [200, 200, 1000, 600]);

% Plot initial and final states
subplot(1,2,1);
plot(x_fm, V_MeV, 'k-', 'LineWidth', 2);
hold on;
plot(x_fm, prob_density(:, 1) * max(V_MeV), 'b-', 'LineWidth', 2);
xlabel('Position (fm)');
ylabel('Energy (MeV)');
title('Initial State');
legend('Potential', 'Probability Density', 'Location', 'best');
grid on;

subplot(1,2,2);
plot(x_fm, V_MeV, 'k-', 'LineWidth', 2);
hold on;
plot(x_fm, prob_density(:, end) * max(V_MeV), 'r-', 'LineWidth', 2);
xlabel('Position (fm)');
ylabel('Energy (MeV)');
title('Final State');
legend('Potential', 'Probability Density', 'Location', 'best');
grid on;

sgtitle(sprintf('Tunneling Animation Summary (E = %.1f MeV)', energy_MeV), ...
        'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'figs/tunneling_summary.png');
saveas(gcf, 'figs/tunneling_summary.fig');

%% Display Results
fprintf('  Animation completed!\n');
fprintf('    Video saved: %s\n', video_file);
fprintf('    Summary plot: figs/tunneling_summary.png\n');
fprintf('    Animation duration: %.1f seconds\n', duration);
fprintf('    Wave packet energy: %.1f MeV\n', energy_MeV);

end
