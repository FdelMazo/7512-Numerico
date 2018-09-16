function integral = integral(lim_inf, lim_sup)
  quad("f",lim_inf,lim_sup)
end

function y = f(x)
  padron1 = 100029; padron2 = 99779;
  P = (padron1 + padron2) / 50;
  alpha = 0.17; beta = 0.41;
  y = ( sin(P*x) + beta * (x^2) ) / (alpha*x);
end