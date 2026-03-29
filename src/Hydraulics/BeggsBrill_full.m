function out = BeggsBrill_full(qL, qG, D, theta, rho_l, rho_g, mu_l, mu_g, sigma)
% BeggsBrill_full  Beggs and Brill correlation with full pattern map
% Inputs:
%   qL, qG - liquid and gas volumetric flow rates [m3/s]
%   D      - pipe inner diameter [m]
%   theta  - inclination [rad] (0 horizontal, pi/2 vertical up)
%   rho_l, rho_g - densities [kg/m3]
%   mu_l, mu_g   - viscosities [Pa.s]
%   sigma - surface tension [N/m]
% Output struct out with fields:
%   FlowPattern, HL, rho_mix, vs_l, vs_g, v_m, Re_m, f, dp_dz

% Guards and defaults
if nargin<9, sigma = 0.02; end
if nargin<8, mu_g = 1e-5; end
if nargin<7, mu_l = 1e-3; end
qL = max(qL,0); qG = max(qG,0);
A = pi*D^2/4;
vs_l = qL / max(A,1e-12);
vs_g = qG / max(A,1e-12);
v_m = vs_l + vs_g;

% superficial velocity ratio and mixture density
jr = vs_l / max(vs_g,1e-12);
rho_m = (rho_l*vs_l + rho_g*vs_g) / max(v_m,1e-12);

% Flow pattern map following Beggs & Brill simplified tables
% Use superficial velocities (m/s) and inclination (deg)
theta_deg = theta * 180/pi;
if theta_deg < -10
    % downward steep: treat as distributed
    pattern = 'Distributed';
elseif theta_deg < 10
    % near horizontal: use jr thresholds
    if jr > 1.2
        pattern = 'Segregated';
    elseif jr > 0.5
        pattern = 'Intermittent';
    else
        pattern = 'Distributed';
    end
else
    % upward incline: shift thresholds
    if jr > 1.5
        pattern = 'Segregated';
    elseif jr > 0.6
        pattern = 'Intermittent';
    else
        pattern = 'Distributed';
    end
end

% Holdup correlations per pattern (empirical)
switch pattern
    case 'Segregated'
        HL = 0.98 * (vs_l ./ (vs_l + vs_g + 1e-12)).^0.5;
    case 'Intermittent'
        HL = 0.85 * (vs_l ./ (vs_l + vs_g + 1e-12)).^0.6;
    otherwise % Distributed
        HL = 0.65 * (vs_l ./ (vs_l + vs_g + 1e-12)).^0.4;
end
HL = min(max(HL,1e-4),0.999);

% Mixture density and Reynolds
rho_mix = HL * rho_l + (1-HL) * rho_g;
mu_mix = HL * mu_l + (1-HL) * mu_g;
Re_m = rho_mix * v_m * D / max(mu_mix,1e-12);

% friction factor: laminar or turbulent (Colebrook approximation via explicit Churchill)
if Re_m < 2300
    f = 64 / max(Re_m,1e-12);
else
    % Churchill explicit approximation for f
    A_ch = (2.457*log(1/((7/Re_m)^0.9 + 0.27*(wellRoughness(D)/D))))^-16;
    B_ch = (37530/Re_m)^16;
    f = 8*((8/Re_m)^12 + 1/(A_ch + B_ch)^(3/2))^(1/12);
    if ~isfinite(f) || f<=0, f = 0.3164 / (Re_m^0.25); end
end

% pressure gradient components
g = 9.80665;
dp_dz_fric = f * rho_mix * v_m^2 / (2*D);
dp_dz_grav = rho_mix * g * sin(theta);
dp_dz_acc = 0;

dp_dz = dp_dz_fric + dp_dz_grav + dp_dz_acc;

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

function r = wellRoughness(~)
% default roughness placeholder (m)
r = 1e-5;
end
