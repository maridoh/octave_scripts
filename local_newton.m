
function [X, Iter] = local_newton(func, x0, kmax, epsilon)
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
    d = H \ (-g);
    x = x + d;
    g = numgradient(fname,{x})';
  endfor
  Iter = k;
endfunction
 
