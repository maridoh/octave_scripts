function [X, Iter] = global_newton(func, x0, kmax, sigma, beta, epsilon, rho, p)
  X = [];
  x = x0;
  fname = func2str(func);
  g = numgradient(fname,{x})';
  for k = 1 : kmax
    X(:,k) = x;
    if (norm(g) <= epsilon)
      break;
    endif
    H = numhessian(fname,{x});
    [d, R] = linsolve(H,-g);
    g = numgradient(fname,{x})';
    if ((R = 0) || (g'*d > -rho*norm(d)^p))
      d = -g;
    endif
    l = 0;
    while (true)
      v1 = func(x + (beta^l)*d);
      v2 = func(x) + sigma*(beta^l)*g'*d;
      if (v1 <= v2)
        break;
      endif
      l = l + 1;
    endwhile
    x = x + (beta^l)*d;
  endfor
  Iter = k;
endfunction