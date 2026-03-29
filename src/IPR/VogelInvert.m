function [qmax, q_at_P] = VogelInvert(Pwf, Pr, qwf)
% VogelInvert  invert Vogel IPR to obtain q_max and evaluate Vogel at pressures
% Inputs:
%   Pwf  - measured bottomhole flowing pressure [Pa]
%   Pr   - reservoir (pseudo) pressure [Pa]
%   qwf  - measured liquid volumetric rate at Pwf (same units as qmax)
% Outputs:
%   qmax   - Vogel q_max (liquid rate at Pwf = 0)
%   q_at_P - function handle @(P) giving Vogel(P,Pr,qmax)

% Safety checks
if Pr <= 0
    error('Pr must be positive.');
end
if Pwf < 0
    error('Pwf must be non-negative.');
end
% Vogel fractional form: q/qmax = 1 - 0.2*(Pwf/Pr) - 0.8*(Pwf/Pr)^2
r = Pwf ./ Pr;
frac = 1 - 0.2*r - 0.8*r.^2;
frac = max(frac, 1e-12); % guard
qmax = qwf ./ frac;
% ensure finite positive
if ~isfinite(qmax) || qmax <= 0
    qmax = max(qwf, 1e-12);
end

% return evaluator
q_at_P = @(P) qmax .* max(0, 1 - 0.2.*(P./Pr) - 0.8.*(P./Pr).^2);
end
