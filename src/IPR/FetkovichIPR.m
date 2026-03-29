function q = FetkovichIPR(t, qi, Di, b)
% FetkovichIPR  Fetkovich type curve approximation
% t  - time since start [days]
% qi - initial rate [m3/day]
% Di - decline rate parameter [1/day]
% b  - decline exponent (0 for exponential)
% returns q(t)
if b==0
    q = qi .* exp(-Di .* t);
else
    q = qi ./ ((1 + b .* Di .* t).^(1./b));
end
end
