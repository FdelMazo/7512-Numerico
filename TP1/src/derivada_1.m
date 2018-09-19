function y = derivada_1(x)
  padron1 = 100029; padron2 = 99779;
  P = (padron1 + padron2) / 50;
  ALPHA = 0.17; BETA = 0.41;
  
  primer_term = (P./(ALPHA.*x)) .* cos(P.*x);
  segundo_term = - ( ( sin(P.*x) ) ./ ( ALPHA .* (x.^2) ) );
  tercer_term = BETA / ALPHA;
  y = primer_term + segundo_term + tercer_term;
end
