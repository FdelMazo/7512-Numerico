function integral = main
  padron1 = 100029; padron2 = 99779;
  global P = (padron1 + padron2) / 50
  global ALPHA = 0.17; global BETA = 0.41;
  global ERR_MAX = 10e-5;
  global LIM_INF = 1; global LIM_SUP = 240; 
  global MUS = 10e-8;
  n = calcular_n(ERR_MAX)
  global n = 10000000;
  integral_d = 6.9458e+04;
  integral_s = 7.0236e+04;
  te = abs(calcular_te(integral_d, integral_s))
  err = te.*MUS
  tic
  n = 10000000
  integral = calcular_area(n)
  toc
  for i = 1:16;
    perturbacion = 1/(10 .^i)
    cp = perturbar(perturbacion,n)
   end
end

function te = calcular_te(d,s)
  global MUS
  te = (d.-s)./(d.*(MUS))
end

function cp = perturbar(perturbacion,n)
  global ALPHA BETA
  ALPHA += perturbacion; BETA += perturbacion;
  valor_perturbado_sup = calcular_area(n)
  ALPHA -= 2 .*perturbacion; BETA -= 2 .*perturbacion;
  valor_perturbado_inf = calcular_area(n)
  ALPHA += perturbacion; BETA += perturbacion;
  cp =  abs((1 .- (valor_perturbado_inf ./ valor_perturbado_sup)) ./ perturbacion);
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
  y = abs(primer_term) + abs(segundo_term) + abs(tercer_term);
end

function y = fderivada2(x)
  global P ALPHA BETA
 
  primer_term =  - (2.*P.*cos(P.*x) ) ./ (ALPHA .* (x.^2) );
  segundo_term = 2* sin(P.*x) ./ (ALPHA .* (x.^3) );
  tercer_term =  - ( ( (P^2)*sin(P.*x) ) ./ (ALPHA .* x) );
  y = abs(primer_term) + abs(segundo_term) + abs(tercer_term);
end

function n = calcular_n(error_maximo_truncamiento)
  global LIM_SUP LIM_INF
  num = - ( (LIM_SUP - LIM_INF)^3 ) * fderivada2(1);
  denom = error_maximo_truncamiento * 12;
  n = sqrt(abs(num/denom));
end

function a = calcular_area(n)
  global LIM_SUP LIM_INF;
  h = ( LIM_SUP - LIM_INF ) ./ n;
  f_inicio = f(LIM_INF) / 2;
  f_fin = f(LIM_SUP) / 2;
  f_i = 0;
  for i = 1:n-1;
    f_i = f_i + f(LIM_INF + i*h) * h;
  end
    a = ( f_inicio + f_fin ) * h + f_i;
end

function g = graficar()
  fplot(@f, [-0.02 0.02])
  fplot(@fderivada, [-0.02 0.02])
  fplot(@fderivada2, [-0.02 0.02])
end

function y = graficar_n()
  x = [1,10,100,1000,10000,100000,1000000,10000000];
  y = [];
  for n = x
    n
    tic;
    integral = calcular_area(n)
    y = [y,toc]; 
    printf("Tiempo = %ds\n\n",y(length(y)))   
  end
  plot(x,y,'o-r')
end