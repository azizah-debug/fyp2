function Z = PR_eos(P, T, comp, x, phase)
% PR_eos  Peng-Robinson cubic solver (robust to comp shape)
if nargin<5, phase='auto'; end
R = 8.314462618; % J/mol/K

% ensure x is row vector
x = reshape(x,1,[]);

% extract component arrays
Tc = reshape([comp.Tc],1,[]);
Pc = reshape([comp.Pc],1,[]);
omega = reshape([comp.omega],1,[]);

nc = numel(Tc);
if numel(x) ~= nc
    error('Length of mole fraction vector x (%d) does not match number of components (%d).', numel(x), nc);
end

% compute ai(T) and bi (row vectors)
m = 0.37464 + 1.54226.*omega - 0.26992.*omega.^2;
alpha = (1 + m.*(1 - sqrt(T./Tc))).^2;
ai = 0.45724 .* R.^2 .* Tc.^2 ./ Pc .* alpha;   % 1 x nc
bi = 0.07780 .* R .* Tc ./ Pc;                 % 1 x nc

% mixing rules (kij = 0) - correct double-sum vectorized
ai_col = ai(:);            % nc x 1
Aij = sqrt(ai_col * ai_col'); % nc x nc
Xij = (x(:) * x(:)');         % nc x nc
a_mix = sum(sum(Xij .* Aij)); % scalar

b_mix = sum(x .* bi);         % scalar

A = a_mix * P / (R^2 * T^2);
B = b_mix * P / (R * T);

% cubic coefficients Z^3 + c2 Z^2 + c1 Z + c0 = 0
c2 = -(1 - B);
c1 = A - 3*B^2 - 2*B;
c0 = -(A*B - B^2 - B^3);
coeffs = [1 c2 c1 c0];
Zroots = roots(coeffs);

% select real roots
Zreal = real(Zroots(abs(imag(Zroots))<1e-8));
if isempty(Zreal)
  Z = real(Zroots(1));
else
  switch lower(phase)
    case 'vapor'
      Z = max(Zreal);
    case 'liquid'
      Z = min(Zreal);
    otherwise
      if numel(Zreal)>1
        Z = max(Zreal);
      else
        Z = Zreal(1);
      end
  end
end
end
