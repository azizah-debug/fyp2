function result = NodalNewton(surfaceP, T, z, comp, qL_surface, qG_surface, well, opts)
% NodalNewton  Newton solver with adaptive acceptance of refined analytic Jacobian
% Usage:
%   result = NodalNewton(surfaceP, T, z, comp, qL_surface, qG_surface, well, opts)
% opts fields (all optional):
%   tolP      - pressure tolerance [Pa] (default 1e2)
%   maxit     - max Newton iterations (default 30)
%   fd_eps    - finite difference step fraction (default 1e-6)
%   verbose   - logical (default false)
%   useAnalyticJacobian - logical (default false)
%   useRefinedAnalyticJacobian - logical (default false)
%   refinedClampPerSeg - per-segment clamp for refined Jacobian sensitivity (default 1e-2)
%   logIterationFile - optional path to CSV file for per-iteration logging
%
% Returns struct result with fields:
%   BHP, converged, iters, res_norm, details

if nargin < 8, opts = struct(); end

% defaults
if isfield(opts,'tolP'), tolP = opts.tolP; else tolP = 1e2; end
if isfield(opts,'maxit'), maxit = opts.maxit; else maxit = 30; end
if isfield(opts,'fd_eps'), fd_eps = opts.fd_eps; else fd_eps = 1e-6; end
if isfield(opts,'verbose'), verbose = opts.verbose; else verbose = false; end
if isfield(opts,'useAnalyticJacobian'), useAnalytic = opts.useAnalyticJacobian; else useAnalytic = false; end
if isfield(opts,'useRefinedAnalyticJacobian'), useRefined = opts.useRefinedAnalyticJacobian; else useRefined = false; end
if isfield(opts,'refinedClampPerSeg'), refinedClampPerSeg = opts.refinedClampPerSeg; else refinedClampPerSeg = 1e-2; end

% prepare Vogel closure for IPR (assumes VogelInvert exists)
[~, qVogel] = VogelInvert(surfaceP, max(surfaceP+1e5,surfaceP*1.1), qL_surface);

% optional logger
if isfield(opts,'logIterationFile') && ~isempty(opts.logIterationFile)
  try
    logger = IterationLogger(opts.logIterationFile);
  catch
    logger = [];
    if verbose
      warning('Could not create IterationLogger at %s. Continuing without logging.', opts.logIterationFile);
    end
  end
else
  logger = [];
end

% residual function closure
resfun = @(P) nodalResidualWrapper(P, surfaceP, T, z, comp, qVogel, qG_surface, well);

% initial guess
P = surfaceP + 5e5;
F = resfun(P);
if ~isfinite(F)
    if ~isempty(logger), logger.close(); end
    error('Initial residual is non-finite at P=%g', P);
end

res_norm = abs(F);
converged = false;

for it = 1:maxit
    if verbose
        fprintf('Newton iter %d: P=%g, res=%g\n', it, P, F);
    end
    if abs(F) < tolP
        converged = true;
        % log final state
        if ~isempty(logger)
            try
                logger.log(it, P, F, 'Done', NaN, NaN);
            catch
            end
        end
        break;
    end

    % Initialize per-iteration logging variables
    jacSource = 'FD';
    alpha_log = NaN;
    dp_log = NaN;

    % Attempt refined analytic Jacobian first (if enabled)
    dFdP_refined = NaN;
    accepted = false;
    if useRefined
        dFdP_refined = refinedAnalyticJacobian_clamped(P, surfaceP, T, z, comp, qVogel, qG_surface, well, refinedClampPerSeg);
        if isfinite(dFdP_refined) && abs(dFdP_refined) > 1e-16
            % form trial Newton step using refined Jacobian
            dp_trial = -F / dFdP_refined;
            % trial update and evaluate residual
            P_trial = max(P + dp_trial, surfaceP + 1e3);
            F_trial = resfun(P_trial);
            if isfinite(F_trial) && abs(F_trial) < abs(F)
                % accept refined Jacobian for this iteration
                dFdP = dFdP_refined;
                jacSource = 'RefinedAnalytic';
                dp = dp_trial;
                Pnew = P_trial;
                Fnew = F_trial;
                accepted = true;
                alpha_log = 1.0;
                dp_log = dp;
            else
                % reject refined Jacobian for this iteration; will use FD below
                dFdP_refined = NaN;
            end
        end
    end

    % If refined not used/accepted, try analytic (simple) if requested
    if ~accepted
        if useAnalytic
            dFdP_analytic = analyticJacobian(P, surfaceP, T, z, comp, qVogel, qG_surface, well);
            if isfinite(dFdP_analytic) && abs(dFdP_analytic) > 1e-16
                % form trial step using analytic Jacobian
                dp_trial = -F / dFdP_analytic;
                P_trial = max(P + dp_trial, surfaceP + 1e3);
                F_trial = resfun(P_trial);
                if isfinite(F_trial) && abs(F_trial) < abs(F)
                    dFdP = dFdP_analytic;
                    jacSource = 'Analytic';
                    dp = dp_trial;
                    Pnew = P_trial;
                    Fnew = F_trial;
                    accepted = true;
                    alpha_log = 1.0;
                    dp_log = dp;
                else
                    dFdP_analytic = NaN;
                end
            end
        end
    end

    % If neither analytic nor refined accepted, use finite-difference Jacobian
    if ~accepted
        h = max(abs(P)*fd_eps, 1e-6);
        Fp = resfun(P + h);
        if ~isfinite(Fp)
            h = h * 0.1;
            Fp = resfun(P + h);
        end
        if ~isfinite(Fp)
            dFdP = (Fp - F) / max(h,1e-12);
        else
            dFdP = (Fp - F) / h;
        end
        jacSource = 'FD';
        % compute step and prepare trial values via line search below
        if ~isfinite(dFdP) || abs(dFdP) < 1e-12
            dp = sign(-F) * max(1e3, 0.1*abs(P));
        else
            dp = -F / dFdP;
        end
    end

    % If we already have Pnew/Fnew from accepted analytic/refined, skip line search
    if exist('Pnew','var') && exist('Fnew','var') && isfinite(Pnew) && isfinite(Fnew)
        % accepted earlier; dp_log already set
    else
        % backtracking line search for dp
        alpha = 1.0; c = 1e-4; maxls = 10; success = false;
        for ls = 1:maxls
            P_try = P + alpha * dp;
            P_try = max(P_try, surfaceP + 1e3);
            F_try = resfun(P_try);
            if ~isfinite(F_try)
                alpha = alpha * 0.5; continue;
            end
            if abs(F_try) <= (1 - c*alpha)*abs(F) || abs(F_try) < tolP
                success = true;
                Pnew = P_try;
                Fnew = F_try;
                alpha_log = alpha;
                dp_log = alpha * dp;
                break;
            else
                alpha = alpha * 0.5;
            end
        end
        if ~success
            % conservative fallback step
            Pnew = P + 0.1 * dp;
            Pnew = max(Pnew, surfaceP + 1e3);
            Fnew = resfun(Pnew);
            alpha_log = 0.1;
            dp_log = 0.1 * dp;
            if ~isfinite(Fnew)
                % log failure and break
                if ~isempty(logger)
                    try
                        logger.log(it, P, F, jacSource, alpha_log, dp_log);
                    catch
                    end
                end
                break;
            end
        end
    end

    if verbose
        fprintf('  jacobian source: %s, step alpha=%g, |F|->%g\n', jacSource, alpha_log, abs(Fnew));
    end

    % update for next iteration
    P = Pnew;
    F = Fnew;
    res_norm = abs(F);

    % log iteration
    if ~isempty(logger)
        try
            logger.log(it, P, F, jacSource, alpha_log, dp_log);
        catch
        end
    end

    % clear temporary vars to avoid accidental reuse
    if exist('dFdP','var'), clear dFdP; end
    if exist('dFdP_refined','var'), clear dFdP_refined; end
    if exist('dFdP_analytic','var'), clear dFdP_analytic; end
    if exist('dp','var'), clear dp; end
    if exist('Pnew','var'), clear Pnew; end
    if exist('Fnew','var'), clear Fnew; end
    if exist('alpha','var'), clear alpha; end
end

% close logger if open
if ~isempty(logger)
  try
    logger.close();
  catch
  end
end

result.BHP = P;
result.converged = converged;
result.iters = it;
result.res_norm = res_norm;
[~, tubeflow, details] = nodalResidualWrapper(P, surfaceP, T, z, comp, qVogel, qG_surface, well);
result.TubingFlow = tubeflow;
result.details = details;
end

%% ------------------------------------------------------------------------
%% Local helper: nodalResidualWrapper (with per-segment diagnostics)
function [res, tubeflow, details] = nodalResidualWrapper(BHP_local, surfaceP, T, z, comp, qVogel, qG_surface, well)
% compute residual p_surface_pred - surfaceP and return details
% Uses Vogel IPR closure and integrates tubing dp using BeggsBrill_full
% Returns per-segment dp_seg and p_seg in details for leak detection

qL_ipr = qVogel(BHP_local);
if ~isfinite(qL_ipr) || qL_ipr < 0, qL_ipr = 0; end
qG_ipr = qG_surface * (BHP_local / max(surfaceP,1e5));
if ~isfinite(qG_ipr) || qG_ipr < 0, qG_ipr = 0; end

Nseg = 20;
dz = well.L / Nseg;
p_local = BHP_local;
dp_total = 0;
beta = 0; out = []; last_props = struct();

% preallocate per-segment arrays
dp_seg = zeros(1,Nseg);
p_seg = zeros(1,Nseg+1);
p_seg(1) = p_local;

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

    out = BeggsBrill_full(qL_loc, qG_loc, well.D, well.theta, rho_l, rho_g, 1e-3, 1e-5, 0.02);
    dpdz = out.dp_dz;
    if ~isfinite(dpdz), dpdz = 1e6; end
    dp = dpdz * dz;
    dp_seg(k) = dp;
    dp_total = dp_total + dp;
    p_local = p_local - dp;
    p_seg(k+1) = p_local;

    if ~isfinite(dp_total) || ~isfinite(p_local)
        res = 1e9;
        tubeflow = [qL_ipr, qG_ipr];
        details.dp_total = dp_total;
        details.p_surface_pred = NaN;
        details.beta_surface = beta;
        details.BeggsLast = out;
        details.dp_seg = dp_seg;
        details.p_seg = p_seg;
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
details.dp_seg = dp_seg;
details.p_seg = p_seg;
end


%% ------------------------------------------------------------------------
% Local helper: analyticJacobian (conservative fallback analytic)
function dFdP = analyticJacobian(BHP_local, surfaceP, T, z, comp, qVogel, qG_surface, well)
% Conservative analytic Jacobian approximation (simpler fallback)
qL_ipr = qVogel(BHP_local);
qG_ipr = qG_surface * (BHP_local / max(surfaceP,1e5));
A = pi*well.D^2/4;
Nseg = 20; dz = well.L / Nseg; p_local = BHP_local;
sum_dpdP = 0;

for k=1:Nseg
    try
        [beta, xL, yV] = isothermalFlash(p_local, T, comp, z);
        propsL = PR_mixProps(p_local, T, comp, xL, 'liquid');
        propsV = PR_mixProps(p_local, T, comp, yV, 'vapor');
        rho_l = propsL.rho; rho_g = propsV.rho;
    catch
        rho_l = 800; rho_g = 50;
    end
    rho_l = max(rho_l,1e-6); rho_g = max(rho_g,1e-6);

    qL_loc = qL_ipr;
    qG_loc = qG_ipr * (p_local / max(surfaceP,1e3));
    vs_l = qL_loc / max(A,1e-12);
    vs_g = qG_loc / max(A,1e-12);
    v_m = vs_l + vs_g;
    rho_mix = (rho_l*vs_l + rho_g*vs_g) / max(v_m,1e-12);

    out = BeggsBrill_full(qL_loc, qG_loc, well.D, well.theta, rho_l, rho_g, 1e-3, 1e-5, 0.02);
    dpdz = max(out.dp_dz,0);

    drho_g_dP = rho_g / max(p_local,1e-6);
    dqGloc_dP = qG_ipr / max(surfaceP,1e-6);
    dvs_g_dP = dqGloc_dP / max(A,1e-12);
    dv_m_dP = dvs_g_dP;

    if v_m <= 0
        d_dpdz_drho = 0;
        d_dpdz_dvm = 0;
    else
        d_dpdz_drho = dpdz / max(rho_mix,1e-12);
        d_dpdz_dvm  = 2 * dpdz / max(v_m,1e-12);
    end

    drho_mix_dP = (rho_g * dvs_g_dP + vs_g * drho_g_dP) / max(v_m,1e-12) - rho_mix * (dv_m_dP / max(v_m,1e-12));
    ddpdz_dP = d_dpdz_drho * drho_mix_dP + d_dpdz_dvm * dv_m_dP;

    sum_dpdP = sum_dpdP + ddpdz_dP * dz;

    dp = dpdz * dz;
    p_local = p_local - dp;
    if ~isfinite(p_local), break; end
end

dFdP = 1 - sum_dpdP;
if ~isfinite(dFdP) || abs(dFdP) < 1e-16
    dFdP = NaN;
end
end

%% ------------------------------------------------------------------------
% Local helper: refinedAnalyticJacobian_clamped
function dFdP = refinedAnalyticJacobian_clamped(BHP_local, surfaceP, T, z, comp, qVogel, qG_surface, well, perSegClamp)
% Refined analytic Jacobian with per-segment clamp parameter
qL_ipr = qVogel(BHP_local);
qG_ipr = qG_surface * (BHP_local / max(surfaceP,1e5));
A = pi*well.D^2/4;
Nseg = 20; dz = well.L / Nseg; p_local = BHP_local;
sum_dpdP = 0;

% conservative viscosity defaults
mu_l_default = 1e-3;
mu_g_default = 1e-5;

for k=1:Nseg
    try
        [beta, xL, yV] = isothermalFlash(p_local, T, comp, z);
        propsL = PR_mixProps(p_local, T, comp, xL, 'liquid');
        propsV = PR_mixProps(p_local, T, comp, yV, 'vapor');
        rho_l = propsL.rho; rho_g = propsV.rho;
    catch
        rho_l = 800; rho_g = 50;
    end
    rho_l = max(rho_l,1e-6); rho_g = max(rho_g,1e-6);

    qL_loc = qL_ipr;
    qG_loc = qG_ipr * (p_local / max(surfaceP,1e3));
    vs_l = qL_loc / max(A,1e-12);
    vs_g = qG_loc / max(A,1e-12);
    v_m = vs_l + vs_g;
    rho_mix = (rho_l*vs_l + rho_g*vs_g) / max(v_m,1e-12);

    out = BeggsBrill_full(qL_loc, qG_loc, well.D, well.theta, rho_l, rho_g, mu_l_default, mu_g_default, 0.02);
    dpdz = max(out.dp_dz,0);

    % drho_g/dP via small FD on PR_mixProps
    dp_fd = max(1e-3 * p_local, 1e2);
    try
        propsV_p = PR_mixProps(p_local + dp_fd, T, comp, yV, 'vapor');
        rho_g_p = propsV_p.rho;
        drho_g_dP = (rho_g_p - rho_g) / dp_fd;
    catch
        drho_g_dP = rho_g / max(p_local,1e-6);
    end

    Re_m = out.Re_m;
    if ~isfinite(Re_m) || Re_m <= 0
        Re_m = max(1, rho_mix * v_m * well.D / mu_l_default);
    end

    if Re_m < 2300
        df_dRe = -64 / (Re_m^2);
    else
        df_dRe = -0.25 * 0.3164 * Re_m^(-1.25);
    end

    mu_mix = max(mu_l_default * 0.5 + mu_g_default * 0.5, 1e-6);
    dqGloc_dP = qG_ipr / max(surfaceP,1e-6);
    dvs_g_dP = dqGloc_dP / max(A,1e-12);
    dv_m_dP = dvs_g_dP;

    if v_m <= 0
        drho_mix_dP = 0;
    else
        drho_mix_dP = (rho_g * dvs_g_dP + vs_g * drho_g_dP) / max(v_m,1e-12) - rho_mix * (dv_m_dP / max(v_m,1e-12));
    end

    dRe_dP = (drho_mix_dP * v_m + rho_mix * dv_m_dP) * well.D / max(mu_mix,1e-12);

    g = 9.80665;
    d_dpdz_drho = (out.f * v_m^2 / (2*well.D)) + g * sin(well.theta);
    d_dpdz_dvm  = out.f * rho_mix * v_m / (well.D);
    d_dpdz_df   = rho_mix * v_m^2 / (2*well.D);

    ddpdz_dP = d_dpdz_drho * drho_mix_dP + d_dpdz_dvm * dv_m_dP + d_dpdz_df * df_dRe * dRe_dP;

    % clamp per-segment sensitivity using provided perSegClamp
    ddpdz_dP = max(min(ddpdz_dP, perSegClamp), -perSegClamp);

    sum_dpdP = sum_dpdP + ddpdz_dP * dz;

    dp = dpdz * dz;
    p_local = p_local - dp;
    if ~isfinite(p_local), break; end
end

dFdP = 1 - sum_dpdP;
if ~isfinite(dFdP) || abs(dFdP) > 1e6 || abs(dFdP) < 1e-12
    dFdP = NaN;
end
end
