function [limiar_Pma]= max_Pmasucc( SNR,lambda, S , p)

syms limiar;

% --- Derivada da probabilidade de Acesso ---%   
D_Pma_succ = (p*(erf((S - limiar + S*lambda)/(2*(S + 2*S*lambda)^(1/2)))/2 - 1/2) + 1)*(p/(2*pi^(1/2)*exp((S - limiar + S*lambda)^2/(4*(S + 2*S*lambda)))*(S + 2*S*lambda)^(1/2)) - (p - 1)/(2*pi^(1/2)*S^(1/2)*exp((S - limiar)^2/(4*S)))) + (p*(p*(erf((S - limiar + S*lambda)/(2*(S + 2*S*lambda)^(1/2)))/2 - 1/2) - (erf((S - limiar)/(2*S^(1/2)))/2 - 1/2)*(p - 1)))/(2*pi^(1/2)*exp((S - limiar + S*lambda)^2/(4*(S + 2*S*lambda)))*(S + 2*S*lambda)^(1/2));

solutions = double(solve(D_Pma_succ , limiar));

[X,Y] = size(solutions);
if( X == 0)
    disp('Sem Solu??o - Pmasucc');
    limiar_Pma = 0;
else
    solutions = abs(solutions);

    % ---Calculo da Pma---%
    Pfa = qfunc((solutions - S)./(sqrt(2*S)));              %% Probabilidade de falso alarme Pf(yn>gama|H0) => simulada
    Pd = qfunc((solutions - S - S * lambda)./(sqrt(2*S + 4*S * lambda)));        %% Probabilidade de detec??o Pf(yn>gama|H1) => simulada
    Pma = (1-p).*(1-Pfa) + p.*(1-Pd);
    Pi = p.*(1-Pd);
    Pma_succ = (1-Pi).* Pma;
    

    [X,Y] = max(Pma_succ);
    limiar_Pma = solutions(Y);
end;


end

