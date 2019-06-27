function [U] = euler(fprima, a, b, u0, h)
  M = (b - a)/h; N = length(u0);
  T = a:h:b;
  U = zeros(M+1,N);
  U(1,:) = u0;
  
  for j = 1:M
      U(j+1,:) = U(j,:) + h*feval(fprima, U(j,:), T(j));
  end
end