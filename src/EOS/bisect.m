function root = bisect(f, a, b, tol, maxit)
if nargin<4, tol=1e-8; end
if nargin<5, maxit=100; end
fa = f(a); fb = f(b);
if fa==0, root=a; return; end
if fb==0, root=b; return; end
if fa*fb > 0
  a = max(0, a-1e-3); b = min(1, b+1e-3);
  fa = f(a); fb = f(b);
  if fa*fb > 0
    root = 0.5*(a+b); return;
  end
end
for it=1:maxit
  c = 0.5*(a+b);
  fc = f(c);
  if abs(fc) < tol || (b-a)/2 < tol
    root = c; return;
  end
  if fa*fc <= 0
    b = c; fb = fc;
  else
    a = c; fa = fc;
  end
end
root = 0.5*(a+b);
end