function mu_simple
  mu=single(1); digitos=single(1); x=single(2);
  while (x>1)
    digitos = digitos+1;
    mu = mu/10;
    x = 1+mu;
  endwhile  
  mu
endfunction