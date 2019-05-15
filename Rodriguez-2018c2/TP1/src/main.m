function integral = main
  padron1 = 100029; padron2 = 99779;
  global P = (padron1 + padron2) / 50;
  global LIM_INF = 1; global LIM_SUP = 240; 
  global ALPHA = 0.17; global BETA = 0.41;
  global ERR_MAX = 10e-5;
  
  n_sin_truncar = calcular_n(ERR_MAX)
  n = 5000
  integral = calcular_area(n)
  cps = calcular_cps(n);
  cpa = cps(1), cpb = cps(2)
  error_redondeo = calcular_err()
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

function n = calcular_n(error_maximo_truncamiento)
  global LIM_SUP LIM_INF
  num = - ( (LIM_SUP - LIM_INF)^3 ) * fderivada2(1);
  denom = error_maximo_truncamiento * 12;
  n = sqrt(abs(num/denom));
end

##### TERMINO DE ESTABILIDAD Y CONDICION DEL PROBLEMA #####

function err = calcular_err()
  mu_single = 1.0000e-08;
  integral_d = 6.9458e+04;
  integral_s = 7.0236e+04;
  te = (integral_d.-integral_s)./(integral_d.*(mu_single));
  te = abs(te)
  err = te.*mu_single;
end

function cps = calcular_cps(n)
  cps_a = []; cps_b = [];
  for i = 1:16;
    perturbacion = 1/(10 .^i);
    cps_a = [cps_a, perturbarA(perturbacion,n)];
    cps_b = [cps_b, perturbarB(perturbacion,n)];
   end
   cps_a
   cps_b
   #plot(1:16,cps_a,'o-r')
   #plot(1:16,cps_b,'o-r')
   cpa = max(cps_a);
   cpb = max(cps_b);
   cps = [cpa,cpb];
end

##### PERTURBACIONES #####

function cp = perturbarA(perturbacion,n)
  global ALPHA
  ALPHA += perturbacion;
  valor_perturbado_sup = calcular_area(n);
  ALPHA -= 2 .*perturbacion;
  valor_perturbado_inf = calcular_area(n);
  ALPHA += perturbacion;
  cp =  abs((1 .- (valor_perturbado_inf ./ valor_perturbado_sup)) ./ perturbacion);
end

function cp = perturbarB(perturbacion,n)
  global BETA
  BETA += perturbacion;
  valor_perturbado_sup = calcular_area(n);
  BETA -= 2 .*perturbacion;
  valor_perturbado_inf = calcular_area(n);
  BETA += perturbacion;
  cp =  abs((1 .- (valor_perturbado_inf ./ valor_perturbado_sup)) ./ perturbacion);
end

##### FUNCION Y SUS DERIVADAS #####

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