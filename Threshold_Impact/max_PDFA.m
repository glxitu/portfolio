function [limiar_PDFA]= max_PDFA( SNR,lambda, S , p)

syms limiar;
   
% --- Derivada da probabilidade ---%   
D_PDFA = (erf((S - limiar + S*lambda)/(2*(S + 2*S*lambda)^(1/2)))/2 + 1/2)/(2*pi^(1/2)*S^(1/2)*exp((S - limiar)^2/(4*S))) + (erf((S - limiar)/(2*S^(1/2)))/2 - 1/2)/(2*pi^(1/2)*exp((S - limiar + S*lambda)^2/(4*(S + 2*S*lambda)))*(S + 2*S*lambda)^(1/2));
   
solutions = double(solve(D_PDFA, limiar));

[X,Y] = size(solutions);
if( X == 0)
    disp('Sem Solu??o - PDFA');
    limiar_PDFA = 0;
else
    solutions = abs(solutions);

    % ---Calculo da Pma---%
    Pfa = qfunc((solutions - S)./(sqrt(2*S)));              %% Probabilidade de falso alarme Pf(yn>gama|H0) => simulada
    Pd = qfunc((solutions - S - S * lambda)./(sqrt(2*S + 4*S * lambda)));        %% Probabilidade de detec??o Pf(yn>gama|H1) => simulada
    PD_PFA = Pd.*(1 - Pfa);

    [X,Y] = max(PD_PFA);
    limiar_PDFA = solutions(Y);
end;

end

