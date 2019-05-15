function valor_taylor = aproximacion_por_taylor(valor, orden_taylor = 2) 
  # Carga el paquete symbolic para poder derivar funciones nativamente  
  # En Linux se instala desde la terminal --> sudo apt-get install octave-symbolic
  # En Windows desde Octave mismo --> pkg install -forge symbolic
  pkg load symbolic;
  
  # Declara a las variables x t como simbolicas, ya que queremos que la funcion
  # se derive genericamente y no para un valor pre-hecho
  syms x t;
  
  
  # Muestra el valor de la funcion evaluada
  valor_real = f(valor)
  
  # Genera la funcion simbolica de taylor
  funcion_simbolica = taylor(@f, orden_taylor);

  # Reemplaza los valores de la funcion simbolica por los recibidos por parametro
  funcion_simbolica = subs(funcion_simbolica, t, 0); # 0 por enunciado

  # De ser de punto flotante, el valor hay que precisarlo (a 2 decimales)
  valor = vpa(valor,2);
  funcion_taylor_evaluada = subs(funcion_simbolica, x, valor);
  
  # Convierte el resultado simbolico a double
  valor_taylor = double(funcion_taylor_evaluada);

  ### Gráficos ###
  intervalo = linspace(0,1.5);
  y1 = f(intervalo);
  y2 = double(subs(funcion_simbolica, x, vpa(intervalo,2)));
  figure;
  plot(intervalo,y1,intervalo,y2);
  title('Comparacion entre Taylor y funcion real');
  legend({'y1 = f(x)','y2 = taylor(x)'},'Location','southwest');
end

function funcion = taylor(f, orden_taylor)
  # Nuevamente x t como simbolicos
  # (Octave no se lleva bien con las variables globales...)
  syms x t;
  
  funcion = f(t);  
  for i = 1:orden_taylor
    # diff, del paquete symbolic, recibe por parametro la funcion a derivar, 
    # de que variable derivarla, y el orden de derivacion
    # Devuelve una función simbolica (la derivada de la f recibida)
    f_derivada = diff(f(t), t, i);
    funcion += (x - t)^i * f_derivada / factorial(i);  
  end
end

function y = f(x)
  y = e.^x.*cos(x); 
end