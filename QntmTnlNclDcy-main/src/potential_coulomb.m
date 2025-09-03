function V = potential_coulomb(params)
%POTENTIAL_COULOMB Generate Coulomb barrier potential for alpha decay
%
%   This function creates a Coulomb barrier potential that models the
%   nuclear potential experienced by an alpha particle during alpha decay.
%   The potential combines a nuclear well with a Coulomb repulsion barrier.
%
%   Inputs:
%       params - Structure containing simulation parameters
%                Required fields: numerical.x, potential.coulomb_radius,
%                                potential.coulomb_charge, potential.barrier_height
%
%   Outputs:
%       V - Potential energy array [J] at each grid point
%           Same size as params.numerical.x
%
%   Potential Form:
%       V(r) = -V0                    for r < R (nuclear well)
%       V(r) = Z1*Z2*e^2/(4πε₀*r)    for r > R (Coulomb barrier)
%       where R is the nuclear radius
%
%   Physical Background:
%       This potential models the alpha decay process where an alpha particle
%       must tunnel through the Coulomb barrier created by the nuclear charge.
%       The Gamow factor describes the tunneling probability.
%
%   Example:
%       params = init_params('barrier_type', 'coulomb');
%       V = potential_coulomb(params);
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

required_fields = {'coulomb_radius', 'coulomb_charge', 'barrier_height'};
for i = 1:length(required_fields)
    if ~isfield(params.potential, required_fields{i})
        error('params.potential must contain %s field', required_fields{i});
    end
end

%% Extract Parameters
x = params.numerical.x;
R = params.potential.coulomb_radius;  % Nuclear radius [m]
Z1 = params.potential.coulomb_charge; % Alpha particle charge [e]
Z2 = 90;  % Typical heavy nucleus charge (e.g., Th-232, U-238)
V0 = params.potential.barrier_height_J;  % Nuclear well depth [J]

%% Physical Constants (Nuclear Units)
% Use nuclear units for better numerical stability and physical interpretation
hbar_c = 197.327;  % MeV·fm (reduced Planck constant times c)
e_sq_over_4pi_eps0 = 1.4397;  % MeV·fm (e²/(4πε₀) in nuclear units)
k_coulomb = Z1 * Z2 * e_sq_over_4pi_eps0;  % Coulomb barrier strength [MeV·fm]

% Convert to SI units for consistency with rest of code
k_coulomb_SI = k_coulomb * 1.602e-19 * 1e6 * 1e-15;  % Convert to J·m

%% Initialize Potential Array
V = zeros(size(x));

%% Calculate Potential
% For 1D simulation, we use |x| as the radial coordinate
r = abs(x);

% Nuclear well region (r < R)
nuclear_region = r < R;
V(nuclear_region) = -V0;

% Coulomb barrier region (r >= R)
coulomb_region = r >= R;
V(coulomb_region) = k_coulomb_SI ./ r(coulomb_region);



%% Validation
if any(isnan(V)) || any(isinf(V))
    error('Potential contains NaN or Inf values');
end

if max(V) <= 0
    warning('Maximum potential is not positive - no barrier present');
end

%% Display Information
fprintf('  Coulomb barrier potential created:\n');
fprintf('    Nuclear radius: %.1f fm\n', R * 1e15);
fprintf('    Nuclear well depth: %.1f MeV\n', V0 / (1.602e-19 * 1e6));
fprintf('    Coulomb barrier height at R: %.1f MeV\n', ...
        k_coulomb / (R * 1e15));
fprintf('    Alpha particle charge: %d e\n', Z1);
fprintf('    Nucleus charge: %d e\n', Z2);

%% Calculate Gamow Factor (for reference)
% The Gamow factor gives the tunneling probability for alpha decay
if isfield(params, 'energy') && isfield(params.energy, 'incident_energy_J')
    E = params.energy.incident_energy_J;
    if E > 0
        % Convert to nuclear units for Gamow factor calculation
        E_MeV = E / (1.602e-19 * 1e6);
        R_fm = R * 1e15;
        
        % Classical turning point in nuclear units
        r_classical_fm = Z1 * Z2 * e_sq_over_4pi_eps0 / E_MeV;
        
        % Gamow factor calculation (nuclear units)
        % G = 2 * integral from R to r_classical of kappa(r) dr
        % where kappa = sqrt(2*m*(V-E))/hbar_c
        alpha_mass_MeV = params.physics.m / (1.602e-19 * 1e6);  % MeV/c²
        
        if r_classical_fm > R_fm
            % Numerical integration of kappa
            r_integration = linspace(R_fm, r_classical_fm, 1000);
            V_integration = Z1 * Z2 * e_sq_over_4pi_eps0 ./ r_integration;
            kappa_integration = sqrt(2 * alpha_mass_MeV * (V_integration - E_MeV)) / hbar_c;
            
            % Only integrate where kappa is real and positive
            valid_idx = kappa_integration > 0;
            if any(valid_idx)
                gamow_factor = exp(-2 * trapz(r_integration(valid_idx), kappa_integration(valid_idx)));
                fprintf('    Gamow factor at %.1f MeV: %.2e\n', E_MeV, gamow_factor);
            end
        end
    end
end

end
