function mu_doble
  mu=1; digitos=1; x=2;
  while (x>1)
    digitos = digitos+1;
    mu = mu/10;
    x = 1+mu;
  endwhile  
  mu
endfunction