function y = main
  padron1 = 100029; padron2 = 99779;
  global P = (padron1 + padron2) / 50;
  global ALPHA = 0.17; global BETA = 0.41;
  global ERR_MAX = 10e-5;
  global LIM_INF = 1; global LIM_SUP = 240; 
  
  n = calcular_n(ERR_MAX)
  y = calcular_area(n)
end

function y = funcion(x)
  global P ALPHA BETA

  y = ( sin(x.*P) + BETA * (x.^2) ) ./ (x.*ALPHA);
end

function y = derivada_1(x)
  global P ALPHA BETA

  primer_term = (P./(ALPHA.*x)) .* cos(P.*x);
  segundo_term = - ( ( sin(P.*x) ) ./ ( ALPHA .* (x.^2) ) );
  tercer_term = BETA / ALPHA;
  y = primer_term + segundo_term + tercer_term;
end

function y = derivada_2(x)
  global P ALPHA BETA
 
  primer_term = - (2.*P.*cos(P.*x) ) ./ (ALPHA .* (x.^2) );
  segundo_term = 2* sin(P.*x) ./ (ALPHA .* (x.^3) );
  tercer_term = - ( ( (P^2)*sin(P.*x) ) ./ (ALPHA .* x) );
  y = primer_term + segundo_term + tercer_term;
end

function n = calcular_n(error_maximo_truncamiento)
  global LIM_SUP LIM_INF
  num = - ( (LIM_SUP - LIM_INF)^3 ) * derivada_2(1);
  denom = error_maximo_truncamiento * 12;
  n = sqrt(abs(num/denom));
end

function a = calcular_area(n)
  global LIM_SUP LIM_INF
  h = ( LIM_SUP - LIM_INF ) / n
  f_inicio = funcion(LIM_INF) / 2
  f_fin = funcion(LIM_SUP) / 2
  f_i = 0
  for i = 1:n-1;
    f_i = f_i + funcion(LIM_INF + i*h) * h
  endfor
  a = ( f_inicio + f_fin ) * h + f_i
end
