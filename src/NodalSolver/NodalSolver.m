function result = NodalSolver(surfaceP, T, z, comp, qL_surface, qG_surface, well)
% NodalSolver  simple nodal steady-state solver
% Inputs:
%   surfaceP      - surface pressure [Pa]
%   T             - temperature [K]
%   z             - overall mole fractions (for EOS/flash)
%   comp          - components struct array (from components.mat)
%   qL_surface    - liquid volumetric rate at surface [m3/s]
%   qG_surface    - gas volumetric rate at surface [m3/s]
%   well          - struct with fields: D [m], L [m], theta [rad], roughness
%
% Output struct result:
%   result.BHP    - computed bottomhole pressure [Pa]
%   result.beta   - vapor fraction at surface (if flash used)
%   result.dp_total - total pressure drop [Pa]
%   result.details - struct with intermediate values

% Basic checks
if nargin < 7
    error('Not enough inputs. Provide surfaceP, T, z, comp, qL_surface, qG_surface, well.');
end

% Initial guess for BHP (surfaceP + hydrostatic)
g = 9.80665;
rho_w = 1000; % crude initial density [kg/m3]
BHP_guess = surfaceP + rho_w * g * well.L * sin(well.theta);

% Define residual function: residual(BHP) = IPR(BHP) - TubingFlow(BHP)
residual = @(BHP) nodalResidual(BHP, surfaceP, T, z, comp, qL_surface, qG_surface, well);

% Solve for BHP using simple bisection between surfaceP and a high upper bound
bl = surfaceP;
bu = surfaceP + 1e8; % 100 MPa upper bound
tol = 1e3; % Pa tolerance
maxit = 60;
for it=1:maxit
    mid = 0.5*(bl+bu);
    rmid = residual(mid);
    rbl = residual(bl);
    if abs(rmid) < tol
        BHP = mid; break;
    end
    if sign(rmid) == sign(rbl)
        bl = mid;
    else
        bu = mid;
    end
    BHP = mid;
end

% Final evaluation
[IPR_val, tubeflow, details] = nodalResidual(BHP, surfaceP, T, z, comp, qL_surface, qG_surface, well);

result.BHP = BHP;
result.IPR = IPR_val;
result.TubingFlow = tubeflow;
result.dp_total = details.dp_total;
result.details = details;
end

%% Helper: nodal residual and flow evaluation
function [res, tubeflow, details] = nodalResidual(BHP, surfaceP, T, z, comp, qL_surf, qG_surf, well)
% For this skeleton:
% - IPR is approximated as a fixed production rate (could be a function of BHP)
% - Tubing pressure drop computed by integrating Beggs & Brill along well

% Simple IPR model: assume surface rates are produced (no reservoir coupling)
IPR_rateL = qL_surf; % m3/s
IPR_rateG = qG_surf; % m3/s

% Integrate tubing pressure drop from BHP to surface by discretizing the well
Nseg = 20;
dz = well.L / Nseg;
p = BHP;
p_local = p;
qL = IPR_rateL;
qG = IPR_rateG;
dp_total = 0;
for k=1:Nseg
    % approximate local pressure and compute local densities via PR_mixProps
    % Use overall composition z and PR_mixProps at local p and T to get rho_g, rho_l approx
    % For simplicity assume phase split at local p using isothermalFlash
    try
        [beta, xL, yV] = isothermalFlash(p_local, T, comp, z);
        propsL = PR_mixProps(p_local, T, comp, xL, 'liquid');
        propsV = PR_mixProps(p_local, T, comp, yV, 'vapor');
        rho_l = propsL.rho;
        rho_g = propsV.rho;
    catch
        % fallback densities
        rho_l = 800;
        rho_g = 50;
    end

    % compute volumetric rates at local conditions (assume incompressible liquid, compressible gas)
    qL_loc = qL; % m3/s (liquid approx)
    qG_loc = qG * (p_local / surfaceP); % crude gas compressibility approx (ideal gas scaling)

    % Beggs & Brill for this segment
    out = BeggsBrill([], qL_loc, qG_loc, well.D, well.theta, rho_l, rho_g, 1e-3, 1e-5, 0.02);
    dp = out.dp_dz * dz; % Pa
    dp_total = dp_total + dp;
    p_local = p_local - dp; % step down towards surface
end

% tubeflow is the pressure at surface predicted by tubing model
p_surface_pred = BHP - dp_total;
tubeflow = [qL, qG]; % rates used
% residual: positive if IPR > TubingFlow (we use simple rate difference)
% Here we return residual as difference between predicted surface pressure and actual surface pressure
res = p_surface_pred - surfaceP;

details.dp_total = dp_total;
details.p_surface_pred = p_surface_pred;
details.beta_surface = beta;
details.BeggsLast = out;
end
