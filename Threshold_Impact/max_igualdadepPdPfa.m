function [limiar_PD_PFA]= max_igualdadepPdPfa( SNR, lambda, S, p)

syms limiar pa;

solutions =[];
L=0;
while (length(solutions) == 0),
    Pfa = qx((limiar - S)/(sqrt(2*S)));              
    Pd = qx((limiar - S - S * lambda)/(sqrt(2*S + 4*S * lambda)));   
    equation1 =  (1-pa)*Pfa - pa*(1-Pd);
    equation1 = subs(equation1,'pa', p);
    solutions = double(solve(equation1,limiar));
    a=20+floor(100*rand);
    digits(a);
    L=L+1;
    disp(L);
end
    
[X,Y] = size(solutions);
if( X == 0 || solutions(1) == -Inf)
     disp('Sem Solu??o - igualdadePdPfa');
     limiar_PD_PFA = 0;

     l = S * lambda;
     limiarmax = qfuncinv(0.1/100) * sqrt(2*S + 4*S*lambda) + S + S*lambda;
     limiarmin = 0;
     limiar = limiarmin:0.01:limiarmax;
     Pfa  = qfunc((limiar - S)./(sqrt(2*S))) ;
     Pd = qfunc((limiar - S - S * lambda)./(sqrt(2*S + 4*S * lambda)));
     C4 = (1-p).*Pfa - p.*(1-Pd);
     loc = find(C4>=0);
     [X,Y] = min(C4(loc)); 
     limiar_PD_PFA = limiar(Y);

else
     % ---Calculo da Pma---%
     Pfa = qfunc((solutions - S)./(sqrt(2*S)));              %% Probabilidade de falso alarme Pf(yn>gama|H0) => simulada
     Pd = qfunc((solutions - S - S * lambda)./(sqrt(2*S + 4*S * lambda)));        %% Probabilidade de detec??o Pf(yn>gama|H1) => simulada
     PD_PFA = p*(1-Pd) - (1-p)*Pfa;

     [X,Y] = max(PD_PFA);
     limiar_PD_PFA = solutions(Y);
end;
 
end

