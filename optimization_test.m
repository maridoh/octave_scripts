1;

#Gradient Rosenbrock function
function y = rosenbrock(x)
   y = (1-x(1))^2 + 100*(x(2) - x(1)^2)^2;
endfunction

#[X, iter] = quasi_newton_method(str2func("rosenbrock"),[-1.2,1]',eye(2),500,10^(-6))
#[X, iter] = global_newton(str2func("rosenbrock"),[-1.2,1]',200,10^(-4),0.5,10^(-6),10^(-8),2.1)
#[X, iter] = local_newton(str2func("rosenbrock"),[-1.2,1]',500,10^(-6))
#[X, iter] = mod_polak_ribiere(str2func("rosenbrock"),[-1.2,1]',500,10^(-4),0.5,10^(-6),0.1,10)

#trigonometric function
function y = trigo(x)
  s = 0;
  for i = 1:4
    s1 = 0;
    for j = 1:4
      s1 = s1 + cos(x(j)) + i*(1 - cos(x(i))) - sin(x(i));
    endfor
    s = s + (4 - s1)^2;
  endfor
  y = s;
endfunction

[X, iter] = quasi_newton_method(str2func("trigo"),[1/4,1/4,1/4,1/4]',eye(4),500,10^(-6))
#[X, iter] = global_newton(str2func("trigo"),[1/4,1/4,1/4,1/4]',200,10^(-4),0.5,10^(-6),10^(-8),2.1)
#[X, iter] = local_newton(str2func("trigo"),[1/4,1/4,1/4,1/4]',500,10^(-6))
#[X, iter] = mod_polak_ribiere(str2func("trigo"),[1/4,1/4,1/4,1/4]',500,10^(-4),0.5,10^(-6),0.1,10)


#Brown function

function y = brown(x)
  y = (x(1) - 10^6)^2 + (x(2) - 2*10^6)^2 + (x(1)*x(2) - 2)^2;
endfunction
#[X, iter] = quasi_newton_method(str2func("brown"),[1;1],eye(2),500,10^(-6));
#[X, iter] = global_newton(str2func("brown"),[1,1]',200,10^(-4),0.5,10^(-6),10^(-8),2.1)
#[X, iter] = local_newton(str2func("brown"),[1,1]',500,10^(-6))
#[X, iter] = mod_polak_ribiere(str2func("brown"),[1,1]',500,10^(-4),0.5,10^(-6),0.1,10)


#Wood function
function y = wood(x)
  y = 100*(x(2) - x(1)^2)^2 + (1-x(1))^2 +90*(x(4) - x(3)^2)^2 + (1 - x(3))^2 + 10*(x(2) + x(4) - 2)^2 + (1/10)*(x(2) - x(4))^2;
endfunction
#[X, iter] = quasi_newton_method(str2func("wood"),[-3;-1;-3;-1],eye(4),50,10^(-6))
#[X, iter] = global_newton(str2func("wood"),[-3;-1;-3;-1],200,10^(-4),0.5,10^(-6),10^(-8),2.1)
#[X, iter] = local_newton(str2func("brown"),[-3;-1;-3;-1]',500,10^(-6))
#[X, iter] = mod_polak_ribiere(str2func("wood"),[-3;-1;-3;-1]',500,10^(-4),0.5,10^(-6),0.1,10)

#
function y = root(x)
  y = sqrt(1+x^2);
endfunction
#[X, iter] = local_newton(str2func("root"),0.5,200,10^(-6))

