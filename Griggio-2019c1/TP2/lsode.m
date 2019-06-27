function [U] = lsode(fprima, a, b, u0, h)
  # Se pone un wrapper sobre la funcion lsode para cambiar su interfaz
  # Así se puede lograr recibir la cantidad de parametros que reciben el resto de los metodos de resolucion numerica
  # Luego, este archivo transforma esos parametros a los que recibe lsode normalmente
  # Finalmente, con ayuda de la funcion 'builtin' de octave, se llama al "verdadero" lsode
  T = a:h:b;
  [U] = builtin('lsode', fprima, u0, T);
end