%% Histograma/Probabilidade da diferen?a entre os finais das tramas dos utilizadores licenciados PU1 e PU2, em "single frame"

clear all;
close all;


%% SU 1
T_PU_ON_1 = 1e5;                    % # tramas PU 1
T_PU_1 = T_PU_ON_1*2;      
mu_on_1 = 5;                        % tempo medio ON PU 1
mu_off_1 = 5;                       % tempo medio OFF PU 1
desvio_1 = 0;                       % desvio PU 1

%% SU 2
T_PU_ON_2 = 1e5;                    % # tramas PU 2
T_PU_2 = T_PU_ON_2*2;      
mu_on_2 = 5;                        % tempo medio ON PU 2
mu_off_2 = 5;                       % tempo medio OFF PU 2
desvio_2 = 0;                       % desvio PU 2

%% Compute Opportunity

gran = 10;
tts = [];
opp = [];
aux = 0;

for (i=1:1:T_PU_ON_1)

%     dev_1 = desvio_1*rand(1,1);
%     tt_on_1 = dev_1+exprnd(mu_on_1,1,1);
%     tt_off_1 = tt_on_1+exprnd(mu_off_1,1,1);
%     
%     dev_2 = desvio_2*rand(1,1);
%     tt_on_2 = dev_2+exprnd(mu_on_2,1,1);  
%     tt_off_2 = tt_on_2+exprnd(mu_off_2,1,1);
    
    dev_1 = desvio_1*rand(1,1);
    tt_on_1 = dev_1+mu_on_1*rand(1,1);
    tt_off_1 = tt_on_1+mu_off_1*rand(1,1);
    
    dev_2 = desvio_2*rand(1,1);
    tt_on_2 = dev_2+mu_on_2*rand(1,1);  
    tt_off_2 = tt_on_2+mu_off_2*rand(1,1);
       
    diff = tt_off_1 - tt_off_2;

    diff = round(diff * gran) / gran;
    
    aux = aux + 1;

    if(isempty(find(tts==diff)))
        tts = [tts diff];
        opp = [opp 1];
    else
        opp(find(tts==diff))=opp(find(tts==diff))+1;
    end

end

opp = opp/(aux);

[tts_sorted, sort_index] = sort(tts);
opp_sorted = opp(sort_index);

plot(tts_sorted, opp_sorted*gran, 'r.-');
grid on
legend('Sim.');

