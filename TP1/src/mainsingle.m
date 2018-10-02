function integral = mainsingle
  padron1 = single(100029); padron2 = single(99779);
  global P = (padron1 + padron2) / single(50)
  global LIM_INF = single(1); global LIM_SUP = single(240); 
  global ALPHA = single(0.17); global BETA = single(0.41);
  global ERR_MAX = single(10e-5);
  
  n_sin_truncar = calcular_n(ERR_MAX)
  n = single(5000)
  integral = calcular_area(n)
end

function a = calcular_area(n)
  global LIM_SUP LIM_INF;
  h = single(( LIM_SUP - LIM_INF ) ./ n);
  f_inicio = single(f(LIM_INF) / 2);
  f_fin = single(f(LIM_SUP) / 2);
  f_i = 0;
  for i = 1:n-1;
    f_i = f_i + f(LIM_INF + i*h) * h;
    f_i = single(f_i);
  end
    a = single(( f_inicio + f_fin ) * h + f_i);
end

function n = calcular_n(error_maximo_truncamiento)
  global LIM_SUP LIM_INF
  
  num = - ( (LIM_SUP - LIM_INF)^3 ) * fderivada2(1);
  num = single(num);
  denom = error_maximo_truncamiento * 12;
  denom = single(denom);
  n = single(sqrt(abs(num/denom)));
end

##### FUNCION Y SUS DERIVADAS #####

function y = f(x)
  x = single(x);
  global P ALPHA BETA

  y = ( single(sin(x.*P)) + BETA * (x.^2) ) ./ (x.*ALPHA);
  y = single(y);
end

function y = fderivada(x)
  x = single(x);
  global P ALPHA BETA

  primer_term = (P./(ALPHA.*x)) .* single(cos(P.*x));
  segundo_term = - ( ( single(sin(P.*x)) ) ./ ( ALPHA .* (x.^2) ) );
  tercer_term = BETA / ALPHA;
  y = abs(single(primer_term)) + abs(single(segundo_term)) + abs(single(tercer_term));
  y = single(y);
end

function y = fderivada2(x)
  x = single(x);
  global P ALPHA BETA
 
  primer_term =  - (2.*P.*single(cos(P.*x)) ) ./ (ALPHA .* (x.^2) );
  segundo_term = 2* single(sin(P.*x)) ./ (ALPHA .* (x.^3) );
  tercer_term =  - ( ( (P^2)*single(sin(P.*x)) ) ./ (ALPHA .* x) );
  y = abs(single(primer_term)) + abs(single(segundo_term)) + abs(single(tercer_term));
  y = single(y);
end