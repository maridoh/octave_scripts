function [X, Iter] = quasi_newton_method-2(fhandle, x0, B0, kmax, epsison)

  X = [];
  x = x0;
  B = B0;
  g = fhandle(x);
  for k = 1 : kmax
    X(:,k) = x;
    if (norm(g) <= epsison)
      break;
    endif
    d = -B*g;
    x_new = x + d;
    s = x_new - x;
    g_new = fhandle(x_new);
    y = g_new - g;
    x = x_new;
    g = g_new;
    B = B + ((s - B*y)*s' + s*(s - B*y)')/(y'*s) - (((s - B*y)'*y)/((y'*s)^2))*(s*s');
  endfor
  Iter = k;

  endfunction
