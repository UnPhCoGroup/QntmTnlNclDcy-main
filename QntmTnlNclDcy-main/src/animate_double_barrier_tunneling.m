function animate_double_barrier_tunneling(params, varargin)
%ANIMATE_DOUBLE_BARRIER_TUNNELING Create double barrier tunneling animation
%
%   This function creates an animated visualization of quantum tunneling through
%   a double barrier potential, similar to the reference visualization.
%   It shows the wave packet approaching, interacting with, and tunneling through
%   the double barrier structure.
%
%   Inputs:
%       params - Structure containing simulation parameters
%       varargin - Optional parameter-value pairs:
%                  'energy' - Energy of the wave packet [MeV] (default: 5)
%                  'width' - Width of the wave packet (default: 2e-15)
%                  'position' - Initial position of wave packet [m] (default: -15e-15)
%                  'duration' - Animation duration in seconds (default: 8)
%                  'fps' - Frames per second (default: 30)
%                  'show_interference' - Show interference effects (default: true)
%
%   Outputs:
%       None (creates and saves animation)
%
%   Method:
%       Creates a Gaussian wave packet and propagates it through the double
%       barrier potential using the time-dependent Schrödinger equation.
%       The animation demonstrates quantum interference effects between
%       the two barriers.
%
%   Example:
%       params = init_params('barrier_type', 'double', 'gap_width', 4e-15);
%       animate_double_barrier_tunneling(params, 'energy', 8, 'duration', 6);
%
%   Author: UniPhi Collective (inspired by reference visualization)
%   Date: July-September 2025

%% Parse Input Arguments
p = inputParser;
addParameter(p, 'energy', 5, @(x) x > 0);  % MeV
addParameter(p, 'width', 2e-15, @(x) x > 0);  % m
addParameter(p, 'position', -15e-15, @isnumeric);  % m
addParameter(p, 'duration', 8, @(x) x > 0);  % seconds
addParameter(p, 'fps', 30, @(x) x > 0);
addParameter(p, 'show_interference', true, @islogical);
parse(p, varargin{:});

energy_MeV = p.Results.energy;
packet_width = p.Results.width;
initial_position = p.Results.position;
duration = p.Results.duration;
fps = p.Results.fps;
show_interference = p.Results.show_interference;

%% Ensure Double Barrier Type
if ~strcmp(params.simulation.barrier_type, 'double')
    warning('Setting barrier_type to "double" for this animation');
    params.simulation.barrier_type = 'double';
end

%% Generate Double Barrier Potential
V = potential_double_barrier(params);

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
fprintf('  Creating double barrier animation (%d frames)...\n', num_frames);

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
video_file = 'figs/double_barrier_tunneling.mp4';
v = VideoWriter(video_file, 'MPEG-4');
v.FrameRate = fps;
open(v);

for frame = 1:num_frames
    clf;
    
    % Main plot: Probability density with potential
    subplot(2,2,1);
    plot(x_fm, V_MeV, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Potential');
    hold on;
    plot(x_fm, prob_density(:, frame) * max(V_MeV), 'b-', 'LineWidth', 2.5, ...
         'DisplayName', '|ψ(x,t)|²');
    xlabel('Position (fm)');
    ylabel('Energy (MeV)');
    title('Double Barrier Tunneling');
    legend('Location', 'best');
    grid on;
    ylim([min(V_MeV), max(V_MeV)]);
    
    % Real part
    subplot(2,2,2);
    plot(x_fm, V_MeV, 'r-', 'LineWidth', 2);
    hold on;
    plot(x_fm, real_part(:, frame) * max(V_MeV) + max(V_MeV)/2, 'g-', 'LineWidth', 2);
    xlabel('Position (fm)');
    ylabel('Energy (MeV)');
    title('Real Part Re(ψ)');
    legend('Potential', 'Re(ψ)', 'Location', 'best');
    grid on;
    ylim([min(V_MeV), max(V_MeV)]);
    
    % Imaginary part
    subplot(2,2,3);
    plot(x_fm, V_MeV, 'r-', 'LineWidth', 2);
    hold on;
    plot(x_fm, imag_part(:, frame) * max(V_MeV) + max(V_MeV)/2, 'm-', 'LineWidth', 2);
    xlabel('Position (fm)');
    ylabel('Energy (MeV)');
    title('Imaginary Part Im(ψ)');
    legend('Potential', 'Im(ψ)', 'Location', 'best');
    grid on;
    ylim([min(V_MeV), max(V_MeV)]);
    
    % Interference analysis (if enabled)
    subplot(2,2,4);
    if show_interference && frame > 1
        % Show probability density evolution
        plot(x_fm, prob_density(:, frame), 'b-', 'LineWidth', 2, ...
             'DisplayName', 'Current');
        hold on;
        plot(x_fm, prob_density(:, max(1, frame-10)), 'b--', 'LineWidth', 1, ...
             'DisplayName', '10 frames ago');
        plot(x_fm, prob_density(:, max(1, frame-20)), 'b:', 'LineWidth', 1, ...
             'DisplayName', '20 frames ago');
        xlabel('Position (fm)');
        ylabel('Probability Density');
        title('Interference Evolution');
        legend('Location', 'best');
        grid on;
    else
        % Show potential structure
        plot(x_fm, V_MeV, 'r-', 'LineWidth', 2);
        xlabel('Position (fm)');
        ylabel('Energy (MeV)');
        title('Double Barrier Structure');
        grid on;
    end
    
    % Add time information
    current_time = time_points(frame) * 1e21;  % Convert to zeptoseconds
    sgtitle(sprintf('Double Barrier Quantum Tunneling - t = %.1f zs (E = %.1f MeV)', ...
                    current_time, energy_MeV), 'FontSize', 14, 'FontWeight', 'bold');
    
    % Write frame to video
    writeVideo(v, getframe(fig));
end

close(v);

%% Create Static Summary Plot
figure('Position', [200, 200, 1200, 600]);

% Plot initial, middle, and final states
subplot(1,3,1);
plot(x_fm, V_MeV, 'r-', 'LineWidth', 2);
hold on;
plot(x_fm, prob_density(:, 1) * max(V_MeV), 'b-', 'LineWidth', 2);
xlabel('Position (fm)');
ylabel('Energy (MeV)');
title('Initial State');
legend('Potential', '|ψ|²', 'Location', 'best');
grid on;

subplot(1,3,2);
plot(x_fm, V_MeV, 'r-', 'LineWidth', 2);
hold on;
plot(x_fm, prob_density(:, round(num_frames/2)) * max(V_MeV), 'b-', 'LineWidth', 2);
xlabel('Position (fm)');
ylabel('Energy (MeV)');
title('Middle State (Interference)');
legend('Potential', '|ψ|²', 'Location', 'best');
grid on;

subplot(1,3,3);
plot(x_fm, V_MeV, 'r-', 'LineWidth', 2);
hold on;
plot(x_fm, prob_density(:, end) * max(V_MeV), 'b-', 'LineWidth', 2);
xlabel('Position (fm)');
ylabel('Energy (MeV)');
title('Final State');
legend('Potential', '|ψ|²', 'Location', 'best');
grid on;

sgtitle(sprintf('Double Barrier Tunneling Summary (E = %.1f MeV)', energy_MeV), ...
        'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'figs/double_barrier_summary.png');
saveas(gcf, 'figs/double_barrier_summary.fig');

%% Display Results
fprintf('  Double barrier animation completed!\n');
fprintf('    Video saved: %s\n', video_file);
fprintf('    Summary plot: figs/double_barrier_summary.png\n');
fprintf('    Animation duration: %.1f seconds\n', duration);
fprintf('    Wave packet energy: %.1f MeV\n', energy_MeV);
fprintf('    Gap width: %.1f fm\n', params.potential.gap_width * 1e15);

end
