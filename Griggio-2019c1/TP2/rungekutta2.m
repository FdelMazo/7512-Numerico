function [U] = rungekutta2(fprima, a, b, u0, h)
  M = (b - a)/h; N = length(u0);
  T = a:h:b;
  U = zeros(M+1,N);
  U(1,:) = u0;
  
  for j = 1:M
    q1 = h*feval(fprima, U(j,:), T(j));
    q2 = h*feval(fprima, U(j,:)+q1, T(j)+h);
    U(j+1,:) = U(j,:) + (q1 + q2)/2;
  end
end