function V = potential_square(params)
%POTENTIAL_SQUARE Generate square barrier potential for quantum tunneling
%
%   This function creates a square barrier potential commonly used in
%   quantum tunneling simulations. The potential consists of a finite
%   square well with a barrier in the center.
%
%   Inputs:
%       params - Structure containing simulation parameters
%                Required fields: numerical.x, potential.barrier_height,
%                                potential.barrier_width, potential.well_width
%
%   Outputs:
%       V - Potential energy array [J] at each grid point
%           Same size as params.numerical.x
%
%   Potential Form:
%       V(x) = 0                    for |x| > well_width/2
%       V(x) = -V0                  for -well_width/2 < x < well_width/2
%              (except barrier region)
%       V(x) = barrier_height       for -barrier_width/2 < x < barrier_width/2
%
%   Example:
%       params = init_params();
%       V = potential_square(params);
%       plot(params.numerical.x * 1e15, V / (1.602e-19 * 1e6));
%
%   Author: Michael (UniPhi Collective)
%   Date: July-September 2025

%% Validate Input
if ~isfield(params, 'numerical') || ~isfield(params.numerical, 'x')
    error('params must contain numerical.x field');
end

if ~isfield(params, 'potential')
    error('params must contain potential field');
end

required_fields = {'barrier_height', 'barrier_width', 'well_width'};
for i = 1:length(required_fields)
    if ~isfield(params.potential, required_fields{i})
        error('params.potential must contain %s field', required_fields{i});
    end
end

%% Extract Parameters
x = params.numerical.x;
barrier_height = params.potential.barrier_height_J;  % Convert to Joules
barrier_width = params.potential.barrier_width;
well_width = params.potential.well_width;

%% Initialize Potential Array
V = zeros(size(x));

%% Define Well Region (negative potential)
well_left = -well_width / 2;
well_right = well_width / 2;
well_region = (x >= well_left) & (x <= well_right);

%% Define Barrier Region (positive potential)
barrier_left = -barrier_width / 2;
barrier_right = barrier_width / 2;
barrier_region = (x >= barrier_left) & (x <= barrier_right);

%% Apply Potential
% Set well potential (negative)
V(well_region) = -barrier_height * 0.1;  % Shallow well, 10% of barrier height

% Override with barrier potential (positive)
V(barrier_region) = barrier_height;

%% Optional: Add smooth transitions (uncomment for smoother potential)
% transition_width = barrier_width * 0.1;  % 10% of barrier width
% 
% % Left transition
% left_transition = (x >= barrier_left - transition_width) & ...
%                   (x < barrier_left);
% V(left_transition) = barrier_height * ...
%     (1 + cos(pi * (x(left_transition) - barrier_left) / transition_width)) / 2;
% 
% % Right transition
% right_transition = (x > barrier_right) & ...
%                    (x <= barrier_right + transition_width);
% V(right_transition) = barrier_height * ...
%     (1 + cos(pi * (x(right_transition) - barrier_right) / transition_width)) / 2;

%% Validation
if any(isnan(V)) || any(isinf(V))
    error('Potential contains NaN or Inf values');
end

if max(V) <= 0
    warning('Maximum potential is not positive - no barrier present');
end

%% Display Information
fprintf('  Square barrier potential created:\n');
fprintf('    Well width: %.1f fm\n', well_width * 1e15);
fprintf('    Barrier width: %.1f fm\n', barrier_width * 1e15);
fprintf('    Barrier height: %.1f MeV\n', barrier_height / (1.602e-19 * 1e6));
fprintf('    Well depth: %.1f MeV\n', abs(V(well_region(1))) / (1.602e-19 * 1e6));

end
