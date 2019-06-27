function [U] = nystrom(fprima, a, b, u0, h)
  M = (b - a)/h; N = length(u0);
  T = a:h:b;
  U = zeros(M+1,N);
  U(1,:) = u0;
  
  # Primera fila por Euler
  U(2,:) = U(1,:) + h*feval(fprima, U(1,:), T(1));

  for j = 2:M+1
    U(j+1,1) = 2*U(j,1) - U(j-1,1) + h^2*feval(fprima, U(j,:), T(j))(2);
    U(j+1,2) = (U(j+1,1) - U(j-1,1))/(2*h);
    U(j+1,3) = 2*U(j,3) - U(j-1,3) + h^2*feval(fprima, U(j,:), T(j))(4);
    U(j+1,4) = (U(j+1,3) - U(j-1,3))/(2*h);
  end
  
  U(M+1,2) = (U(M+1,1) - U(M,1))/(2*h);
  U(M+1,4) = (U(M+1,3) - U(M,3))/(2*h);
end