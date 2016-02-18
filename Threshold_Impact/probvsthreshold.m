function [Pma , Pma_succ, PD_PFA, PDFA, PMAD] = probvsthreshold(SNR, lambda, S , p)

  
%--- Limites ---%
l = S * lambda;
limiarmax = qfuncinv(0.1/100) * sqrt(2*S + 4*S*lambda) + S + S*lambda;
limiarmin = l - l*0.6;

limiar = limiarmin:1:limiarmax;
    
% --- Probabilidades de detec??o e falso alarme
    Pma.Pfa_b = qfunc((limiar - S)./(sqrt(2*S))) ;
    Pma.Pd_b= qfunc((limiar - S - S * lambda)./(sqrt(2*S + 4*S * lambda)));
    Pma.p = ones(size(limiar));
    Pma.p = Pma.p.*(1 - p);
    
% --- Probabilidade Pma ---%
    Pma.SNR = SNR;
    Pma.P_band = (1-p).*(1-Pma.Pfa_b) + p.*(1-Pma.Pd_b); % Probabilidade de Acesso ao Meio
    Pma.limiar = limiar;
% --------------------------%
  
%--- Probabilidade Pma_sucesso ---%
    Pma_succ.Pi_b = p.*(1-Pma.Pd_b);
    Pma_succ.P_band = (1-Pma_succ.Pi_b).*Pma.P_band;
    Pma_succ.limiar = limiar;
    
    [Y,X]= max(Pma_succ.P_band);
    Pma_succ.Pmaxnum = Y;
    Pma_succ.limiarmaxnum = limiar(X);
% -----------------------------%
    
%--- Probabilidade P = p*(1-Pd) - (1-p)*Pfa ---%
    PD_PFA.P_band = p.*(1-Pma.Pd_b) - (1-p).*Pma.Pfa_b;
    PD_PFA.limiar = limiar;
% -----------------------------%
  
%--- Probabilidade P = Pd*(1-Pfa) ---%
    PDFA.P_band = Pma.Pd_b.*(1-Pma.Pfa_b);
    PDFA.limiar = limiar;
    
    [Y,X]= max(PDFA.P_band);
    PDFA.Pmaxnum = Y;
    PDFA.limiarmaxnum = limiar(X);
% -----------------------------%
    
%--- Probabilidade  PMAD = Pd*Pma ---%
    PMAD.P_band = Pma.Pd_b.*Pma.P_band;
    PMAD.limiar = limiar;
    
    [Y,X]= max(PMAD.P_band);
    PMAD.Pmaxnum = Y;
    PMAD.limiarmaxnum = limiar(X);
% -----------------------------%
    
end

