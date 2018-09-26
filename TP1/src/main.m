function integral = main
  padron1 = 100029; padron2 = 99779;
  global P = (padron1 + padron2) / 50
  global ALPHA = 0.17; global BETA = 0.41;
  global ERR_MAX = 10e-5;
  global LIM_INF = 1; global LIM_SUP = 240; 
  n = calcular_n(ERR_MAX)
  integral = calcular_todas_las_areas(n);
end

function y = f(x)
  global P ALPHA BETA

  y = ( sin(x.*P) + BETA * (x.^2) ) ./ (x.*ALPHA);
end

function y = fderivada(x)
  global P ALPHA BETA

  primer_term = (P./(ALPHA.*x)) .* cos(P.*x);
  segundo_term = - ( ( sin(P.*x) ) ./ ( ALPHA .* (x.^2) ) );
  tercer_term = BETA / ALPHA;
  y = primer_term + segundo_term + tercer_term;
end

function y = fderivada2(x)
  global P ALPHA BETA
 
  primer_term = - (2.*P.*cos(P.*x) ) ./ (ALPHA .* (x.^2) );
  segundo_term = 2* sin(P.*x) ./ (ALPHA .* (x.^3) );
  tercer_term = - ( ( (P^2)*sin(P.*x) ) ./ (ALPHA .* x) );
  y = primer_term + segundo_term + tercer_term;
end

function n = calcular_n(error_maximo_truncamiento)
  global LIM_SUP LIM_INF
  num = - ( (LIM_SUP - LIM_INF)^3 ) * fderivada2(1);
  denom = error_maximo_truncamiento * 12;
  n = sqrt(abs(num/denom));
end

function a = calcular_todas_las_areas(n)
  global LIM_INF LIM_SUP;
  resultado = 0;
  base = ( LIM_SUP - LIM_INF ) / n;
  ini = LIM_INF;
  fin = LIM_INF + base;
  for i = 0:n;
    resultado += calcular_area_trapecio(base,ini,fin);
    ini+=base;
    fin+=base;
  end
  a = resultado;
end

function a_trap = calcular_area_trapecio(base,ini,fin)
    if (f(ini)<f(fin)); altura_rectangulo = (f(ini)); else; altura_rectangulo = (f(fin)); end; 
    if (f(ini)>f(fin)); altura_triangulo = (f(ini)); else; altura_triangulo = (f(fin)); end;
    area_rectangulo = base * altura_rectangulo; 
    area_triangulo = (base * altura_triangulo) / 2;
    a_trap = area_rectangulo + area_triangulo;
end

