function out = BeggsBrill(q_m, q_l, q_g, D, theta, rho_l, rho_g, mu_l, mu_g, sigma)
% BeggsBrill  Beggs & Brill multiphase correlation (compact)
% Inputs:
%   q_m  - total mass flow rate [kg/s] (optional, can pass [] and use q_l+q_g)
%   q_l  - liquid volumetric flow rate [m3/s]
%   q_g  - gas volumetric flow rate [m3/s]
%   D    - pipe inner diameter [m]
%   theta- inclination angle [rad] (0 horizontal, pi/2 vertical up)
%   rho_l, rho_g - liquid and gas densities [kg/m3]
%   mu_l, mu_g   - liquid and gas viscosities [Pa.s]
%   sigma - surface tension [N/m]
%
% Outputs struct out with fields:
%   out.FlowPattern - string
%   out.HL          - liquid holdup (fraction)
%   out.rho_mix     - mixture density [kg/m3]
%   out.v_m         - mixture superficial velocity [m/s]
%   out.dp_dz       - pressure gradient [Pa/m] (positive = drop)
%
% Notes: This is a compact engineering implementation of Beggs & Brill.
%       It uses standard correlations for flow pattern, holdup and friction.

% Ensure row scalars
D = double(D); theta = double(theta);

% Basic derived quantities
A = pi*D^2/4;
vs_l = q_l / A;   % superficial liquid velocity [m/s]
vs_g = q_g / A;   % superficial gas velocity [m/s]
v_m = vs_l + vs_g;

% Mixture properties
rho_m = (rho_l*vs_l + rho_g*vs_g) / max(v_m,1e-12);

% Dimensionless groups
g = 9.80665;
Re_l = rho_l * vs_l * D ./ max(mu_l,1e-12);
Re_g = rho_g * vs_g * D ./ max(mu_g,1e-12);
We = rho_g * vs_g^2 * D / max(sigma,1e-12);

% Flow pattern map (simplified Beggs & Brill)
% Use superficial velocity ratio and inclination to classify
jr = vs_l ./ max(vs_g,1e-12);
if vs_g < 0.1 && vs_l > 0.01
    pattern = 'Segregated';
elseif jr > 1.2
    pattern = 'Segregated';
elseif jr > 0.5
    pattern = 'Intermittent';
else
    pattern = 'Distributed';
end

% Liquid holdup correlations (simplified)
switch pattern
    case 'Segregated'
        HL = 0.98 * (vs_l ./ (vs_l + vs_g + 1e-12)).^0.5;
    case 'Intermittent'
        HL = 0.8 * (vs_l ./ (vs_l + vs_g + 1e-12)).^0.6;
    otherwise % Distributed
        HL = 0.6 * (vs_l ./ (vs_l + vs_g + 1e-12)).^0.4;
end
HL = min(max(HL, 1e-4), 0.999);

% Mixture density using holdup
rho_mix = HL * rho_l + (1-HL) * rho_g;

% Friction factor (mixture approach) - use Blasius for turbulent, laminar otherwise
Re_m = rho_mix * v_m * D / ( (mu_l*HL + mu_g*(1-HL)) + 1e-12 );
if Re_m < 2000
    f = 64 / max(Re_m,1e-12);
else
    f = 0.3164 / (Re_m^0.25);
end

% Pressure gradient components
dp_dz_fric = f * rho_mix * v_m^2 / (2*D);         % friction [Pa/m]
dp_dz_grav = rho_mix * g * sin(theta);           % gravity (positive drop)
% acceleration term (approx)
dp_dz_acc = 0; % neglected in this compact implementation

dp_dz = dp_dz_fric + dp_dz_grav + dp_dz_acc;

% Package outputs
out.FlowPattern = pattern;
out.HL = HL;
out.rho_mix = rho_mix;
out.vs_l = vs_l;
out.vs_g = vs_g;
out.v_m = v_m;
out.Re_m = Re_m;
out.f = f;
out.dp_dz = dp_dz;
end
