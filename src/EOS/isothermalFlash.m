function [beta, x, y] = isothermalFlash(P, T, comp, z, tol, maxit)
% isothermalFlash  isothermal flash using Rachford-Rice with Wilson initial K
% Ensures consistent row-vector shapes to avoid implicit expansion issues.

if nargin<5, tol=1e-8; end
if nargin<6, maxit=100; end

% ensure z is row and normalised
z = reshape(z,1,[]);
if abs(sum(z)-1) > 1e-12
    z = z / sum(z);
end

% extract arrays as row vectors
Tc = reshape([comp.Tc],1,[]);
Pc = reshape([comp.Pc],1,[]);
omega = reshape([comp.omega],1,[]);
nc = numel(Tc);
if numel(z) ~= nc
    error('Length of overall composition z (%d) does not match number of components (%d).', numel(z), nc);
end

% initial K via Wilson (row vector)
K = (Pc ./ P) .* exp(5.37*(1+omega) .* (1 - Tc./T));
K(~isfinite(K) | K<=0) = 1e-6;
K = reshape(K,1,[]); % ensure row

% Rachford-Rice function (works with row K)
fRR = @(b, Kvec) sum( z .* (Kvec - 1) ./ (1 + b*(Kvec - 1)) );

% initial beta via bisection with current K
beta = bisect(@(b) fRR(b,K), 0, 1, tol, 200);

% compute phase compositions (all row vectors)
den = (1 + beta*(K - 1));
x = z ./ den;
y = K .* x;

% normalise (safeguard)
x = x ./ sum(x);
y = y ./ sum(y);

% iterate with fugacity-based K update (successive substitution)
for it=1:maxit
  propsL = PR_mixProps(P, T, comp, x, 'liquid');
  propsV = PR_mixProps(P, T, comp, y, 'vapor');
  phiL = reshape(propsL.fugacityCoeffs,1,[]);
  phiV = reshape(propsV.fugacityCoeffs,1,[]);
  phiV(phiV<=0) = 1e-12;
  Knew = (phiL ./ phiV);
  Knew = max(Knew, 1e-12);
  Knew = reshape(Knew,1,[]);
  % solve RR for new beta
  beta_new = bisect(@(b) fRR(b,Knew), 0, 1, tol, 200);
  den_new = (1 + beta_new*(Knew - 1));
  x_new = z ./ den_new;
  y_new = Knew .* x_new;
  x_new = x_new ./ sum(x_new);
  y_new = y_new ./ sum(y_new);
  if max(abs(beta_new - beta)) < 1e-7 && max(abs(x_new - x)) < 1e-7
    beta = beta_new; x = x_new; y = y_new; break;
  end
  beta = beta_new; x = x_new; y = y_new;
end
end
