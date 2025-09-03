function V = potential_double_barrier(params)
%POTENTIAL_DOUBLE_BARRIER Generate double square barrier potential for quantum tunneling
%
%   This function creates a double square barrier potential that demonstrates
%   quantum interference effects, similar to the reference visualization.
%   The potential consists of two separate barriers with a gap between them.
%
%   Inputs:
%       params - Structure containing simulation parameters
%                Required fields: numerical.x, potential.barrier_height,
%                                potential.barrier_width, potential.well_width
%                Optional fields: potential.gap_width (default: 2 * barrier_width)
%
%   Outputs:
%       V - Potential energy array [J] at each grid point
%           Same size as params.numerical.x
%
%   Potential Form:
%       V(x) = 0                    for |x| > well_width/2
%       V(x) = -V0                  for -well_width/2 < x < well_width/2
%              (except barrier regions)
%       V(x) = barrier_height       for barrier regions
%
%   Double Barrier Structure:
%       [Left Barrier] [Gap] [Right Barrier]
%       Each barrier has width = barrier_width
%       Gap has width = gap_width (default: 2 * barrier_width)
%
%   Example:
%       params = init_params('barrier_type', 'double');
%       V = potential_double_barrier(params);
%       plot(params.numerical.x * 1e15, V / (1.602e-19 * 1e6));
%
%   Author: UniPhi Collective (inspired by reference visualization)
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

% Gap width between barriers (default: 2 * barrier_width)
if isfield(params.potential, 'gap_width')
    gap_width = params.potential.gap_width;
else
    gap_width = 2 * barrier_width;  % Default gap
end

%% Initialize Potential Array
V = zeros(size(x));

%% Define Well Region (negative potential)
well_left = -well_width / 2;
well_right = well_width / 2;
well_region = (x >= well_left) & (x <= well_right);

%% Define Double Barrier Structure
% Left barrier
left_barrier_center = -gap_width / 2 - barrier_width / 2;
left_barrier_left = left_barrier_center - barrier_width / 2;
left_barrier_right = left_barrier_center + barrier_width / 2;
left_barrier_region = (x >= left_barrier_left) & (x <= left_barrier_right);

% Right barrier
right_barrier_center = gap_width / 2 + barrier_width / 2;
right_barrier_left = right_barrier_center - barrier_width / 2;
right_barrier_right = right_barrier_center + barrier_width / 2;
right_barrier_region = (x >= right_barrier_left) & (x <= right_barrier_right);

%% Apply Potential
% Set well potential (negative)
V(well_region) = -barrier_height * 0.1;  % Shallow well, 10% of barrier height

% Override with barrier potentials (positive)
V(left_barrier_region) = barrier_height;
V(right_barrier_region) = barrier_height;



%% Validation
if any(isnan(V)) || any(isinf(V))
    error('Potential contains NaN or Inf values');
end

if max(V) <= 0
    warning('Maximum potential is not positive - no barrier present');
end

%% Display Information
fprintf('  Double barrier potential created:\n');
fprintf('    Well width: %.1f fm\n', well_width * 1e15);
fprintf('    Barrier width: %.1f fm each\n', barrier_width * 1e15);
fprintf('    Gap width: %.1f fm\n', gap_width * 1e15);
fprintf('    Barrier height: %.1f MeV\n', barrier_height / (1.602e-19 * 1e6));
fprintf('    Well depth: %.1f MeV\n', abs(V(well_region(1))) / (1.602e-19 * 1e6));
fprintf('    Total barrier separation: %.1f fm\n', gap_width * 1e15);

end
