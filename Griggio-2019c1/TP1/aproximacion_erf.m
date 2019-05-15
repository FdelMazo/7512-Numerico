function aproximacion_erf(x, precision = 5) 
  output_precision(precision)
  valor_erf = erf(x)
  aproximacion_trapecio = trapecio(x)
  aproximacion_simpson = simpson(x) 
  
  intervalo = linspace(0,4);
  y1 = erf(intervalo);
  y2 = trapecio(intervalo);
  y3 = simpson(intervalo);
  figure
  plot(intervalo,y1,intervalo,y2,intervalo,y3)
  title('Métodos de Aproximación')
  legend({'y1 = erf(x)','y2 = trapecio(x)', 'y3= simpson(x)'},'Location','southwest')
end

function aproximacion = simpson(b)
  a = 0;
  m = (a+b) / 2;
  aproximacion = ((b - a) / 6) .* (f(a) + 4 * f(m) + f(b));
end   

function aproximacion = trapecio(b)
  a = 0;
  aproximacion = (b - a) .* (f(a) + f(b)) / 2;
end   

function y = f(t)
    y = ( 2 / sqrt(pi) ) * e.^ -(t.^2);
end