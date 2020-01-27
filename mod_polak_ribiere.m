function [X, Iter] = mod_polak_ribiere(func, x0, kmax, sigma, beta, epsilon, delta1, delta2)
  X = [];
  x = x0;
  fname = func2str(func);
  g = numgradient(fname,{x})';
  d = -g;
  for k = 1 : kmax
    X(:,k) = x;
    if (norm(g) <= epsilon)
      break;
    endif
    rho = abs(g'*d)/norm(d);
    l = 0;
    while (true)
      x_new = x + (beta^l)*d;
      g_new = numgradient(fname,{x_new})';
      beta_pr = (g_new'*(g_new  - g)) / norm(g)^2;
      d_new = -g_new + beta_pr*d;
      v1 = func(x_new);
      v2 = func(x) + sigma*(beta^(2*l))*norm(d)^2;
      v3 = g_new'*d_new;
      norm_g_new = norm(g_new)^2;
      if ((v1 <= v2 ) && (-delta2*norm_g_new <= v3) && (v3 <= -delta1*norm_g_new))
        break;
      endif
      l = l + 1;
    endwhile
    x = x_new;
    g = g_new;
    d = d_new;
  endfor
  Iter = k;
endfunction