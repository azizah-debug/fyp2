function result = NodalSolverNewton(surfaceP, T, z, comp, qL_surface, qG_surface, well)
% NodalSolverNewton  nodal solver using secant iteration and Vogel inversion
% Inputs:
%   surfaceP, T, z, comp - as before
%   qL_surface, qG_surface - measured surface volumetric rates [m3/s]
%   well - struct with D, L, theta, roughness
% Output:
%   result struct with BHP, TubingFlow, details

% compute Vogel q_max from measured operating point (use surfaceP as Pwf)
Pr_est = max(surfaceP + 1e5, surfaceP*1.1); % conservative Pr estimate if not provided
try
  [qmax, qVogel] = VogelInvert(surfaceP, Pr_est, qL_surface);
catch
  % fallback if VogelInvert fails
  qmax = max(qL_surface, 1e-9);
  qVogel = @(P) max(0, qmax .* (1 - 0.2.*(P./Pr_est) - 0.8.*(P./Pr_est).^2));
end

% secant initial guesses
p0 = surfaceP + 1e5;
p1 = surfaceP + 5e5;
res0 = nodalResidual(p0);
res1 = nodalResidual(p1);

tolP = 1e3; maxit = 60;
BHP = p1;
for it=1:maxit
    if ~isfinite(res0) || ~isfinite(res1)
        warning('NodalSolverNewton:nonfiniteResidual','Non-finite residual encountered; aborting.');
        break;
    end
    if abs(res1-res0) < 1e-12
        p_new = 0.5*(p0+p1);
    else
        p_new = p1 - res1*(p1-p0)/(res1-res0);
    end
    p_new = max(p_new, surfaceP + 1e3);
    p_new = min(p_new, surfaceP + 1e8);
    res_new = nodalResidual(p_new);
    if ~isfinite(res_new)
        warning('NodalSolverNewton:resNaN','Residual became non-finite at p=%g; stopping.', p_new);
        break;
    end
    if abs(res_new) < tolP
        BHP = p_new; break;
    end
    p0 = p1; res0 = res1;
    p1 = p_new; res1 = res_new;
    BHP = p_new;
end

[~, tubeflow, details] = nodalResidual(BHP);
result.BHP = BHP;
result.TubingFlow = tubeflow;
result.details = details;

% Nested helper uses qVogel closure
function [res, tubeflow, details] = nodalResidual(BHP_local)
    % compute IPR rates using Vogel with Pr_est
    qL_ipr = qVogel(BHP_local);
    if ~isfinite(qL_ipr) || qL_ipr < 0, qL_ipr = 0; end
    qG_ipr = qG_surface * (BHP_local / max(surfaceP,1e5));
    if ~isfinite(qG_ipr) || qG_ipr < 0, qG_ipr = 0; end

    % integrate tubing pressure drop
    Nseg = 20; dz = well.L / Nseg; p_local = BHP_local; dp_total = 0;
    beta = 0; out = []; last_props = struct();
    for k=1:Nseg
        try
            [beta, xL, yV] = isothermalFlash(p_local, T, comp, z);
            propsL = PR_mixProps(p_local, T, comp, xL, 'liquid');
            propsV = PR_mixProps(p_local, T, comp, yV, 'vapor');
            rho_l = propsL.rho; rho_g = propsV.rho;
        catch
            rho_l = 800; rho_g = 50;
        end
        if ~isfinite(rho_l) || rho_l <= 0, rho_l = 800; end
        if ~isfinite(rho_g) || rho_g <= 0, rho_g = 50; end

        qL_loc = qL_ipr;
        qG_loc = qG_ipr * (p_local / max(surfaceP,1e3));
        qL_loc = max(qL_loc, 0); qG_loc = max(qG_loc, 0);

        out = BeggsBrill([], qL_loc, qG_loc, well.D, well.theta, rho_l, rho_g, 1e-3, 1e-5, 0.02);
        dpdz = out.dp_dz;
        if ~isfinite(dpdz), dpdz = 1e6; end
        dp = dpdz * dz;
        dp_total = dp_total + dp;
        p_local = p_local - dp;
        if ~isfinite(dp_total) || ~isfinite(p_local)
            res = 1e9;
            tubeflow = [qL_ipr, qG_ipr];
            details.dp_total = dp_total;
            details.p_surface_pred = NaN;
            details.beta_surface = beta;
            details.BeggsLast = out;
            return;
        end
    end
    p_surface_pred = BHP_local - dp_total;
    tubeflow = [qL_ipr, qG_ipr];
    res = p_surface_pred - surfaceP;
    details.dp_total = dp_total;
    details.p_surface_pred = p_surface_pred;
    details.beta_surface = beta;
    details.BeggsLast = out;
end
end
