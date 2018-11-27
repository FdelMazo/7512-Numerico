function test = main
  padron = 99779;
  diametro_externo = 244.48 *10^-3;
  global DENSIDAD = 7850;
  global PI = 3.14159265359;
  global CALOR_ESPECIFICO = 480;
  global LARGO_TUBO = 12;
  global LARGO_HORNO = 50;
  global NBOL = 50;
  global temperatura1 = round( 473 / 10000 * (padron - 90000) + 773 );
  global temperatura2 = round( 473 / 10000 * (padron - 90000) + 773 );
  global HC = 20;
  global SIGMA = 5.6703*10^-8;
  global EPSILON = 0.85;
  global t0 = 293;
  global S = diametro_externo * LARGO_TUBO * PI
  cadencia = round( -10  / 10000 * (padron - 90000) + 35  );
  diametro_interno = 0.01384;
  velocidad0 = LARGO_HORNO / (NBOL * cadencia);
  termino1 = DENSIDAD * PI * diametro_externo * diametro_interno * LARGO_TUBO;
  termino2 = 1 - diametro_interno/diametro_externo;
  masa = termino1 * termino2;
  
  tiempo = 0:cadencia:LARGO_HORNO/velocidad0
  
  exactos = metodo_exacto(tiempo, masa)
  euler = metodo_euler(velocidad0, cadencia, masa)
  rk = conveccion_rk(velocidad0, cadencia, masa);
  
  plotear_temps(tiempo, euler, rk, exactos);
  plotear_errores(tiempo, euler, rk, exactos);
endfunction

function void = plotear_errores(t, euler, rk, exactos)
  error_euler = calcular_error(t, exactos, euler);
  error_runge = calcular_error(t, exactos, rk);
  plot(t, error_euler)
  hold on
  plot(t, error_runge)
  hold off
endfunction

function error = calcular_error(t, exactos, otro)
  e = [];
  err = 0;
  for i = 1:columns(exactos);
    e = [e err];
    err += abs(exactos(i) - otro(i));
  endfor
  error = e;
endfunction

function void = plotear_temps(t, euler, rk, exactas)
  t = t ./ 60;
  plot(t, exactas .- 273);
  title("T(t)")
  xlabel("t(m)")
  ylabel("T(C°)")
  hold on
  plot(t, euler .- 273)
  plot(t, rk .-273)
  hold off
endfunction
  
 function temperaturas = metodo_euler(velocidad0, cadencia, masa)
  global LARGO_HORNO
  temp = 293;
  fin = LARGO_HORNO/velocidad0;
  v = [];
  for t = cadencia:cadencia:fin;
    v = [v temp];
    temp += cadencia*conveccion_diferencial(temp, masa);
  endfor
  temperaturas = [v temp];
endfunction

function cc = conveccion_diferencial(temp, masa)
  global temperatura1 HC CALOR_ESPECIFICO S
  termino1 = temp - temperatura1;
  termino2 = -masa * CALOR_ESPECIFICO;
  cc = HC * S * termino1 / termino2;
endfunction

function exactas = metodo_exacto(t, masa)
  for tiempo = t;
    valor = conveccion_exacta(t,masa);
  endfor
  exactas = valor;
endfunction

function cc = conveccion_exacta(t, masa)
  global temperatura1 HC t0 CALOR_ESPECIFICO S
  termino1 = -HC*S*t;
  termino2 = masa * CALOR_ESPECIFICO;
  exponente = termino1/termino2;
  a = exp(termino1/termino2);
  cc = temperatura1 + (t0 - temperatura1)*exp(exponente);
endfunction

function rk = conveccion_rk(velocidad0, cadencia, masa)
global LARGO_HORNO
  temp = 293;
  fin = LARGO_HORNO/velocidad0;
  v = [];
  for t = cadencia:cadencia:fin;
    v = [v temp];
    temp += (cadencia/6)*rk4(temp, masa, cadencia);
  endfor
  rk = [v temp];
endfunction

function k = rk4(temp, masa, h)
  k1 = conveccion_diferencial(temp, masa);
  k2 = conveccion_diferencial(temp + h/2, masa);
  k3 = conveccion_diferencial(temp + h/2, masa);
  k4 = conveccion_diferencial(temp + h, masa);
  k = k1 + 2*k2 + 2*k3 + k4;
endfunction
