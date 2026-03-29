function props = PR_mixProps(P, T, comp, x, phase)
% PR_mixProps  compute Z, Vm, rho, MWmix, fugacityCoeffs, h_molar (approx)
if nargin<5, phase='auto'; end
R = 8.314462618; % J/mol/K

% ensure x is row vector and normalised
x = reshape(x,1,[]);
if abs(sum(x)-1) > 1e-12
    x = x / sum(x);
end

% extract arrays
Tc = reshape([comp.Tc],1,[]);
Pc = reshape([comp.Pc],1,[]);
omega = reshape([comp.omega],1,[]);
Mw = reshape([comp.Mw],1,[]);
nc = numel(Tc);
if numel(x) ~= nc
    error('Length of mole fraction vector x (%d) does not match number of components (%d).', numel(x), nc);
end

% compressibility factor
Z = PR_eos(P, T, comp, x, phase);

% molar volume and mixture MW
Vm = Z * R * T / P; % m3/mol
Mw_gmol = sum(x .* Mw); % g/mol
MWmix = Mw_gmol / 1000; % kg/mol

% density
rho = MWmix / Vm; % kg/m3

% compute ai, bi and mixture a,b for fugacity coefficient calculation
m = 0.37464 + 1.54226.*omega - 0.26992.*omega.^2;
alpha = (1 + m.*(1 - sqrt(T./Tc))).^2;
ai = 0.45724 .* R.^2 .* Tc.^2 ./ Pc .* alpha; % 1 x nc
bi = 0.07780 .* R .* Tc ./ Pc;               % 1 x nc

% correct vectorized a_mix and b_mix using double-sum
ai_col = ai(:);
Aij = sqrt(ai_col * ai_col');   % nc x nc
Xij = (x(:) * x(:)');           % nc x nc
a_mix = sum(sum(Xij .* Aij));
b_mix = sum(x .* bi);

A = a_mix * P / (R^2 * T^2);
B = b_mix * P / (R * T);

% compute fugacity coefficients (PR residual form)
phi = zeros(1,nc);
for i=1:nc
  a_i_mix = sum( x .* sqrt(ai(i) .* ai) ); % scalar
  term1 = bi(i)/b_mix * (Z - 1) - log(max(Z - B,1e-12));
  logarg = (Z + (1+sqrt(2))*B) ./ (Z + (1-sqrt(2))*B);
  term2 = -A/(2*sqrt(2)*B) * (2*a_i_mix/a_mix - bi(i)/b_mix) * log(max(logarg,1e-12));
  lnphi = term1 + term2;
  phi(i) = exp(lnphi);
end

% approximate ideal gas molar enthalpy using simple cp placeholder
cp_comp = 20 + 0.01*Mw; % J/mol/K (placeholder)
h_ideal = sum(x .* cp_comp) * (T - 298.15); % J/mol relative to 298.15K

% enthalpy departure (approx)
hd = R*T*(Z - 1) - a_mix/(2*sqrt(2)*b_mix) * log(max((Z + (1+sqrt(2))*B) / (Z + (1-sqrt(2))*B),1e-12));

h_molar = h_ideal + hd;

props.Z = Z;
props.Vm = Vm;
props.rho = rho;
props.MWmix = MWmix;
props.fugacityCoeffs = phi;
props.h_molar = h_molar;
end
