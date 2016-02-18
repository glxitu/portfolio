function [ limiar_Pma, limiar_Pmasucc, limiar_equalpPdPfa, limiar_equalPdPfa, limiar_Pdfa, limiar_Pmad, limiar_TNP] = threshold_computation( SNR, lambda, S, p, interfFU5)

    % --- Calculo do threshold para Maximiza??o da Probabilidade de Acesso
    %[limiar_Pma]= max_Pma( SNR , lambda, S , p);  
    limiar_Pma = -1;
    
    % --- Calculo do threshold para Maximiza??o da Probabilidade de Acesso
    % com sucesso
    %[limiar_Pmasucc]= max_Pmasucc(SNR, lambda, S , p);  %  necessita de max Num?rico
    limiar_Pmasucc = -1; 
    
    % --- Calculo do threshold para Maximizar p*(1-Pd) = (1-p)*(1-Pfa)
    [limiar_equalpPdPfa]= max_igualdadepPdPfa( SNR, lambda, S , p);
    %limiar_equalpPdPfa = -1;
    
    % --- Calculo do threshold para Maximizar (1-Pd) = (1-Pfa)
    %[limiar_equalPdPfa]= max_igualdadePdPfa( SNR, lambda, S);
    limiar_equalPdPfa = -1;
    
    % --- Calculo do threshold para Maximizar PDFA = Pd*(1-Pfa)
    %[limiar_Pdfa]= max_PDFA( SNR, lambda, S , p); % necessita de max Num?rico
    limiar_Pdfa = -1;
        
    % --- Calculo do threshold para Maximizar PMAD = Pd*Pma
    %[ limiar_Pmad ] = max_PMAD( SNR, lambda, S , p); % necessita de max Num?rico
    limiar_Pmad = -1;
    
    % Teorema Neyman-Pearson
    %[ limiar_TNP ] = Neyman_Pearson( SNR,lambda, S , p, interfFU5);  
    limiar_TNP = -1;
    
     if(limiar_Pmasucc == 0 || limiar_Pmad==0 || limiar_Pdfa==0)
         disp('Erro : Solve');
         [Pma, Pma_succ, PD_PFA, PDFA, PMAD] = probvsthreshold(SNR, lambda, S , p);
         if(limiar_Pmasucc==0)
             limiar_Pmasucc = Pma_succ.limiarmaxnum
         end;
         if(limiar_Pmad == 0)
            limiar_Pmad = PMAD.limiarmaxnum
         end;
         if(limiar_Pdfa==0)
            limiar_Pdfa = PDFA.limiarmaxnum
         end;
     end;    
end

