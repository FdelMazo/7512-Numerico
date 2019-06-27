function [U] = rungekutta4(fprima, a, b, u0, h)
  M = (b - a)/h; N = length(u0);
  T = a:h:b;
  U = zeros(M+1,N);
  U(1,:) = u0;
  
  for j = 1:M  
    q1 = h*feval(fprima, U(j,:), T(j));
    q2 = h*feval(fprima, U(j,:) + q1/2, T(j) + h/2);
    q3 = h*feval(fprima, U(j,:) + q2/2, T(j) + h/2);
    q4 = h*feval(fprima, U(j,:) + q3, T(j) + h/2);    
    U(j+1,:) = U(j,:) + (q1 + 2*q2 + 2*q3 + q4)/6;
  end
end