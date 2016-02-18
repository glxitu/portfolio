function [limiar]= Neyman_Pearson(SNR,lambda, S , p, interfFU5)

Pd = 1 - (interfFU5/(100*p));
limiar = qfuncinv(Pd) * sqrt(2*S + 4*S*lambda) + S + S*lambda;  %% Limiar de pot?ncia (gama)
   
end
