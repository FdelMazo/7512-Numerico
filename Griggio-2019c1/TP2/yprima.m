function [y] = yprima(y0, t)
  # y0 = [x1, x1', x2, x2'] = [x1, v1, x2, v2]
  # y' = [x1', x1'', x2', x2''] = [v1, v1', v2, v2']
  constantes; # Carga del archivo consantes.m mu y eta
  d1 = norm( [y0(1)+mu, y0(3)] );
  d2 = norm( [y0(1)-eta, y0(3)] );
  
  y(1) = y0(2);
  y(2) = 2*y0(4) + y0(1) - eta*((y0(1)+mu)/d1^3) - mu*((y0(1)-eta)/d2^3);
  y(3) = y0(4);
  y(4) = -2*y0(2) + y0(3) - eta*(y0(3)/d1^3) - mu*(y0(3)/d2^3);
end
