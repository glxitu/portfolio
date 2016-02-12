%% Histograma/Probabilidade da diferen?a entre os finais das tramas dos utilizadores licenciados PU1 e PU2, em sequ?ncia

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

gran=10;
%% 

% tt_on_1 = exprnd(mu_on_1,T_PU_ON_1,1);
% tt_off_1 = exprnd(mu_off_1,T_PU_ON_1,1);

tt_on_1 = mu_on_1*rand(T_PU_ON_1,1);
tt_off_1 = mu_off_1*rand(T_PU_ON_1,1);

tt_1=[];
for(i=1:1:T_PU_ON_1)
    tt_1 = [tt_1 tt_on_1(i) tt_off_1(i)];
end
tt_1=cumsum(tt_1);

% tt_on_2 = exprnd(mu_on_2,T_PU_ON_2,1);
% tt_off_2 = exprnd(mu_off_2,T_PU_ON_2,1);

tt_on_2 = mu_on_2*rand(T_PU_ON_2,1);
tt_off_2 = mu_off_2*rand(T_PU_ON_2,1);

tt_2=[];
for(i=1:1:T_PU_ON_2)
    tt_2 = [tt_2 tt_on_2(i) tt_off_2(i)];
end
tt_2=cumsum(tt_2);

gran = 10;
diff = [];
quant_diff = [];
total = 0;

for (i=1:2:2*T_PU_ON_1)
    
    diff_aux = tt_1(i)-tt_2(i);
    
    diff_aux = round(diff_aux * gran) / gran;
    
    if(isempty(find(diff==diff_aux)))
        diff = [diff diff_aux];
        quant_diff = [quant_diff 1];
    else
        quant_diff(find(diff==diff_aux))=quant_diff(find(diff==diff_aux))+1;
    end
end

quant_diff=quant_diff/T_PU_ON_1;

[diff_sorted, sort_index] = sort(diff);
quant_diff_sorted = quant_diff(sort_index);

plot(diff_sorted, quant_diff_sorted, 'r.-');
grid on
legend('Sim.');

