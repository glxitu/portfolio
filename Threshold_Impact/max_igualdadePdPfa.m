function [limiar_PD_PFA]= max_igualdadePdPfa( SNR, lambda, S)

syms limiar;

Pfa = qx((limiar - S)/(sqrt(2*S)));              
Pd = qx((limiar - S - S * lambda)/(sqrt(2*S + 4*S * lambda)));   
equation1 =  Pfa - (1-Pd);
solutions = double(solve(equation1,limiar));
 
[X,Y] = size(solutions);
if( X == 0)
    disp('Sem Solu??o - igualdadePdPfa');
    limiar_PD_PFA = 0;
else
     % ---Calculo da Pma---%
     Pfa = qfunc((solutions - S)./(sqrt(2*S)));              %% Probabilidade de falso alarme Pf(yn>gama|H0) => simulada
     Pd = qfunc((solutions - S - S * lambda)./(sqrt(2*S + 4*S * lambda)));        %% Probabilidade de detec??o Pf(yn>gama|H1) => simulada
     PD_PFA = (1-Pd) - Pfa;

     [X,Y] = max(PD_PFA);
     limiar_PD_PFA = solutions(Y);
end;

end

