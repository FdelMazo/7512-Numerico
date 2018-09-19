function y = funcion(x)
  padron1 = 100029; padron2 = 99779;
  P = (padron1 + padron2) / 50;
  ALPHA = 0.17; BETA = 0.41;

  y = ( sin(x.*P) + BETA * (x.^2) ) ./ (x.*ALPHA);
end