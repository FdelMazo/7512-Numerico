function y = derivada_2(x)
  padron1 = 100029; padron2 = 99779;
  P = (padron1 + padron2) / 50;
  ALPHA = 0.17; BETA = 0.41;
  
  primer_term = - (2.*P.*cos(P.*x) ) ./ (ALPHA .* (x.^2) );
  segundo_term = 2* sin(P.*x) ./ (ALPHA .* (x.^3) );
  tercer_term = - ( ( (P^2)*sin(P.*x) ) ./ (ALPHA .* x) );
  y = primer_term + segundo_term + tercer_term;
end
