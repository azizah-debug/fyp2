function q = VogelIPR(Pwf, Pr, q_max)
% VogelIPR  Vogel inflow performance relationship for solution-gas drive oil wells
% Pwf  - bottomhole flowing pressure [Pa]
% Pr   - reservoir (or pseudo) pressure [Pa]
% q_max - maximum liquid rate at Pwf=0 (same units as q)
% Returns q (liquid volumetric rate)

% Convert to relative pressure (Vogel uses fraction)
if Pr <= 0
    error('Pr must be positive');
end
Prf = Pwf ./ Pr;
% Vogel correlation (fractional form)
q = q_max .* (1 - 0.2.*Prf - 0.8.*Prf.^2);
% enforce non-negative
q(q<0) = 0;
end
