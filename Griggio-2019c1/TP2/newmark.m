function [U] = newmark(fprima, a, b, u0, h)
  M = (b - a)/h; N = length(u0);
  T = a:h:b;
  U = zeros(M+1,N);
  U(1,:) = u0;
 
  for j = 1:M
    U(j+1,1) = U(j,1) + h*U(j,2) + h^2 /2 * feval(fprima, U(j,:), T(j))(2);
    U(j+1,3) = U(j,3) + h*U(j,4) + h^2 /2 * feval(fprima, U(j,:), T(j))(4);
    U(j+1,2) = U(j,2) + h/2 * (feval(fprima, U(j,:), T(j))(2) + feval(fprima, [U(j+1,1), U(j,:)(2:end)], T(j+1))(2));
    U(j+1,4) = U(j,4) + h/2 * (feval(fprima, U(j,:), T(j))(4) + feval(fprima, [U(j,:)(1:2), U(j+1,3), U(j,4)], T(j+1))(4));   
  end

end