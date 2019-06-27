function y = main()
  valoresiniciales; #Carga del archivo valoresiniciales.m los valores iniciales del problema
  y0 = [x1 v1 x2 v2];
  op = input("Elija el método para resolver las ecuaciones diferenciales: \n\
    1.lsode \n\
    2.Euler \n\
    3.Runge-Kutta (orden 2)\n\
    4.Runge-Kutta (orden 4)\n\
    5.Nyström \n\
    6.Newmark \n\n\
Opción: ");
  
  orden = 2; # Por defecto se llama a rungekutta tiene orden 2
  if     (op == 1) metodo = 'lsode'; f = @lsode;
  elseif (op == 2) metodo = 'Euler';f = @euler;
  elseif (op == 3) metodo = 'Runge-Kutta (orden 2)'; f = @rungekutta2;
  elseif (op == 4) metodo = 'Runge-Kutta (orden 4)'; f = @rungekutta4;
  elseif (op == 5) metodo = 'Nystrom'; f = @nystrom;
  elseif (op == 6) metodo = 'Newmark';f = @newmark;
  end
  
  [y] = f('yprima',t0, t1, y0, h, orden);  

  figure
  plot(y(:,1), y(:,3))
  title(['x2 en funcion de x1 -- Metodo ' metodo])
  
end