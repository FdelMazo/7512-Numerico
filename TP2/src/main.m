function test = main
  padron = 99779;
  diametro_externo = 244.48 *10^-3;
  global DENSIDAD = 7850;
  global PI = 3.14159265359;
  global CALOR_ESPECIFICO = 480;
  global LARGO_TUBO = 12;
  global LARGO_HORNO = 50;
  global NBOL = 50;
  temperatura1 = round( 200 / 10000 * (padron - 90000) + 500 );
  temperatura2 = round( 200 / 10000 * (padron - 90000) + 500 );
  temperatura1 = temperatura1 + 273;
  temperatura2 = temperatura2 + 273;
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
  fin = LARGO_HORNO/velocidad0;
  
  tiempo = 0:cadencia:LARGO_HORNO/velocidad0;
  
  # Calculo de ED
  exactos = metodo_exacto(tiempo, masa, temperatura1);
  euler = metodo_euler(velocidad0, cadencia, masa, temperatura1);
  rk = conveccion_rk(velocidad0, cadencia, masa, temperatura1);
  
  # Grafico de errores y metodos Euler / RK4
  plotear_temps(tiempo, euler, rk, exactos);
  plotear_errores(tiempo, euler, rk, exactos);
  
  # Comparacion radiacion y conveccion
  r_rk4 = conveccion_radacion_rk(fin, cadencia, masa, temperatura1, t0);
  r_rk4 = r_rk4.-273;
  plotear_comparacion(tiempo, r_rk4, rk) 
  
  # Soaking
  sk = calcular_indice_soaking(r_rk4);
  tiempo_soaking = (fin - tiempo(sk)) / 60;
  temp_soaking = calcular_temp_soaking(r_rk4, sk);
  printf("El tiempo de soaking para T1 = %d °C y T2 = %d °C es:\n%d minutos\n\n", temperatura1, temperatura2, tiempo_soaking);
  printf("Su temperatura de soaking es:\n%d °C\n\n", temp_soaking);
  
  # Tanteo
  t1 = 780;
  t2 = 680;
  rk = rk_dividido(fin, cadencia, masa, t1+273, t2+273, t0);
  rk = rk.-273;
  sk = calcular_indice_soaking(rk);
  temp_soaking = calcular_temp_soaking(rk, sk);
  tiempo_soaking = (fin - tiempo(sk)) /60;
  printf("Se tantearon valores hasta llegar a T1 = %d y T2 = %d para conseguir un tiempo de soaking de 10 minutos\n", t1, t2);
  
  # Calculo automatico
  jacobiano = inv([0.75, 0.25 ; 0.25 , 0.75 ]);
  tsk_obj1 = round( 100 / 10000 * (padron - 90000) + 550 );
  tsk_obj2 = round( 100 / 10000 * (padron - 90000) + 600 );
  tisk_obj = 10;
  tempsA = obtener_temps(fin, cadencia, masa, jacobiano, [667.45;667.45], tisk_obj, 667.45, tiempo);
  tempsB = obtener_temps(fin, cadencia, masa, jacobiano, [tsk_obj1;tsk_obj1], tisk_obj, tsk_obj1, tiempo);
  tempsC = obtener_temps(fin, cadencia, masa, jacobiano, [tsk_obj2;tsk_obj2], tisk_obj, tsk_obj2, tiempo);
  rkA = rk_dividido(fin, cadencia, masa, tempsA(1)+273, tempsA(2)+273, t0);
  rkB = rk_dividido(fin, cadencia, masa, tempsB(1)+273, tempsB(2)+273, t0);
  rkC = rk_dividido(fin, cadencia, masa, tempsC(1)+273, tempsC(2)+273, t0);
  printf("\nCalculos automaticos (10 minutos de soaking para todas):\n\n")
  printf("Las temperaturas de inicio para una temperatura de soaking de %d °C son:\n", temp_soaking);
  printf("T1 = %d °C\nT2 = %d °C\n\n", tempsA(1), tempsA(2));  
  printf("Las temperaturas de inicio para una temperatura de soaking de %d °C son:\n", tsk_obj1);
  printf("T1 = %d °C\nT2 = %d °C\n\n", tempsB(1), tempsB(2));
  printf("Las temperaturas de inicio para una temperatura de soaking de %d °C son:\n", tsk_obj2);
  printf("T1 = %d °C\nT2 = %d °C\n\n", tempsC(1), tempsC(2));
  
  plot(tiempo./60, rkA.-273)
  plot(tiempo./60, rkB.-273)
  plot(tiempo./60, rkC.-273)
  title("T(t)")
  xlabel("t(m)")
  ylabel("T(C°)")
endfunction

function temps = obtener_temps(fin, cadencia, masa, jacobiano, v0, tisk_obj, tsk_obj, tiempo)
  err_individual = 1;
  i = 0;
  while err_individual > 0.5*10^-2;
    v1 = punto_fijo(fin, cadencia, masa, jacobiano, v0, tisk_obj, tsk_obj, tiempo);
    res_err = [(v1(1) - v0(1)) / v1(1); (v1(2) - v0(2)) / v1(2)];
    err_individual = norm(res_err, "inf");
    v0 = v1;
    i+=1;
    if (i == 200)
      break
    endif
  endwhile
  i = i;
  temps = v1;
endfunction

function temps = punto_fijo(fin, h, masa, jacobiano, temperaturas, tiemsk, tempsk, tiempo)
  rk = rk_dividido(fin, h, masa, temperaturas(1) + 273, temperaturas(2) + 273);
  rk = rk.-273;
  sk = calcular_indice_soaking(rk);
  temp_soaking = calcular_temp_soaking(rk, sk);
  tiempo_soaking = (fin - tiempo(sk)) / 60;
  v = [tiempo_soaking - tiemsk; temp_soaking - tempsk];
  temps = temperaturas - jacobiano * v;
endfunction

function rk = rk_dividido(fin, h, masa, t1, t2)
  v = conveccion_radacion_rk(fin/2, h, masa, t1, 273);
  len = columns(v);
  v2 = conveccion_radacion_rk(fin/2, h, masa, t2, v(len));
  for i = 2:len
    v = [v v2(i)];  
  endfor
  rk = v;
endfunction

function void = plotear_comparacion(tiempo, r_rk4, rk)
  plot(tiempo./60, r_rk4)
  hold on
  title("T(t)")
  xlabel("t(m)")
  ylabel("T(C°)")
  legend("Radiacion + Conveccion")
  plot(tiempo./60, rk.-273)
  hold off
endfunction
  
function indice_soaking = calcular_indice_soaking(r_rk4)
  len = columns(r_rk4);
  sk = 0;
  for i = 1:columns(r_rk4);
        if ((r_rk4(len) - r_rk4(i)) > 10)
          sk = i+1;
        endif
  endfor
  indice_soaking = sk;
endfunction

function temp_soaking = calcular_temp_soaking(r_rk4, sk)
  temp_soaking = 0;
  for i = sk:columns(r_rk4)
    temp_soaking += r_rk4(i);
  endfor
  temp_soaking = temp_soaking / (columns(r_rk4) - sk + 1);
  
endfunction

function void = plotear_errores(t, euler, rk, exactos)
  error_euler = calcular_error(exactos, euler);
  error_runge = calcular_error(exactos, rk);
  stem(t ./ 60, log10(error_euler))
  title("e(t)")
  xlabel("t(m)")
  ylabel("error")
  legend("Error Euler")
  hold on
  stem(t ./ 60, log10(error_runge))
  hold off
endfunction

function error = calcular_error(exactos, otro)
  e = [];
  exactos = exactos .- 273;
  otro = otro .- 273;
  err = 0;
  #e = [e err]
  for i = 1:columns(exactos);
    err = (exactos(i) - otro(i)) / exactos(i);
    e = [e err];
  endfor
  error = e;
endfunction

function void = plotear_temps(t, euler, rk, exactas)
  t = t ./ 60;
  ex = plot(t, exactas .- 273);
  title("T(t)")
  xlabel("t(m)")
  ylabel("T(C°)")
  legend(ex, "Valor exacto")
  hold on
  eu = plot(t, euler .- 273);
  rk = plot(t, rk .-273);
  hold off
endfunction
  
 function temperaturas = metodo_euler(velocidad0, cadencia, masa, tinf)
  global LARGO_HORNO
  temp = 293;
  fin = LARGO_HORNO/velocidad0;
  v = [];
  for t = cadencia:cadencia:fin;
    v = [v temp];
    temp += cadencia*conveccion_diferencial(temp, masa, tinf);
  endfor
  temperaturas = [v temp];
endfunction

function cc = conveccion_diferencial(temp, masa, tinf)
  global HC CALOR_ESPECIFICO S
  termino1 = temp - tinf;
  termino2 = -masa * CALOR_ESPECIFICO;
  cc = HC * S * termino1 / termino2;
endfunction

function exactas = metodo_exacto(t, masa, tinf)
  for tiempo = t;
    valor = conveccion_exacta(t,masa, tinf);
  endfor
  exactas = valor;
endfunction

function cc = conveccion_exacta(t, masa, tinf)
  global HC t0 CALOR_ESPECIFICO S
  termino1 = -HC*S*t;
  termino2 = masa * CALOR_ESPECIFICO;
  exponente = termino1/termino2;
  a = exp(termino1/termino2);
  cc = tinf + (t0 - tinf)*exp(exponente);
endfunction

function rk = conveccion_rk(velocidad0, cadencia, masa, tinf)
global LARGO_HORNO
  temp = 293;
  fin = LARGO_HORNO/velocidad0;
  v = [];
  for t = cadencia:cadencia:fin;
    v = [v temp];
    temp += (cadencia/6)*rk4(temp, masa, cadencia, tinf);
  endfor
  rk = [v temp];
endfunction

function k = rk4(temp, masa, h, tinf)
  k1 = conveccion_diferencial(temp, masa, tinf);
  k2 = conveccion_diferencial(temp + h/2, masa, tinf);
  k3 = conveccion_diferencial(temp + h/2, masa, tinf);
  k4 = conveccion_diferencial(temp + h, masa, tinf);
  k = k1 + 2*k2 + 2*k3 + k4;
endfunction

function rk = conveccion_radacion_rk(fin, cadencia, masa, tinf, t0)
  v = [];
  temp = t0;
  for t = cadencia:cadencia:fin;
    v = [v temp];
    temp += (cadencia/6)*r_rk4(temp, masa, cadencia, tinf);
  endfor
  rk = [v temp];
endfunction

function k = r_rk4(temp, masa, h, tinf)
  k1 = conveccion_radacion(temp, masa, tinf);
  k2 = conveccion_radacion(temp + h/2, masa, tinf);
  k3 = conveccion_radacion(temp + h/2, masa, tinf);
  k4 = conveccion_radacion(temp + h, masa, tinf);
  k = k1 + 2*k2 + 2*k3 + k4;
endfunction

function cc = conveccion_radacion(temp, masa, tinf)
  global SIGMA EPSILON CALOR_ESPECIFICO S
  divisor = - masa * CALOR_ESPECIFICO;
  numerador = SIGMA * EPSILON * S * (temp^4 - tinf^4);
  radiacion = numerador / divisor;
  cc = conveccion_diferencial(temp, masa, tinf) + radiacion;
endfunction
