%% Quantifica??o da oportunidade de transmiss?o dispon?vel considerando a utiliza??o do canal por parte de 2 utilizadores licenciados, PU1 e PU2

clear all;
close all;

%% SU 1
T_PU_ON_1 = 1e3;                    % # tramas PU 1
T_PU_1 = T_PU_ON_1*2;      
mu_on_1 = 5;                        % tempo medio ON PU 1
mu_off_1 = 5;                       % tempo medio OFF PU 1
desvio_1 = 0;                       % desvio PU 1

%% SU 2
T_PU_ON_2 = 1e3;                    % # tramas PU 2
T_PU_2 = T_PU_ON_2*2;      
mu_on_2 = 5;                        % tempo medio ON PU 2
mu_off_2 = 5;                       % tempo medio OFF PU 2
desvio_2 = 0;                       % desvio PU 2

gran=10;
%% 

tt_on_1 = exprnd(mu_on_1,T_PU_ON_1,1);
tt_off_1 = exprnd(mu_off_1,T_PU_ON_1,1);

% tt_on_1 = mu_on_1*rand(T_PU_ON_1,1);
% tt_off_1 = mu_off_1*rand(T_PU_ON_1,1);

tt_1=[];
for(i=1:1:T_PU_ON_1)
    tt_1 = [tt_1 ones(1,round(tt_on_1(i)*gran)) zeros(1,round(tt_off_1(i)*gran))];
end

tt_on_2 = exprnd(mu_on_2,T_PU_ON_2,1);
tt_off_2 = exprnd(mu_off_2,T_PU_ON_2,1);

% tt_on_2 = mu_on_2*rand(T_PU_ON_2,1);
% tt_off_2 = mu_off_2*rand(T_PU_ON_2,1);

tt_2=[];
for(i=1:1:T_PU_ON_2)
    tt_2 = [tt_2 ones(1,round(tt_on_2(i)*gran)) zeros(1,round(tt_off_2(i)*gran))];
end

excesso_trama = 0;
if length(tt_1)>length(tt_2)
    excesso_trama=length(tt_1)-length(tt_2);
    tt_1 = tt_1(1:length(tt_2));
else
    excesso_trama=length(tt_2)-length(tt_1);
    tt_2 = tt_2(1:length(tt_1));
end

opp = not(tt_1 | tt_2);
opp_perc = mean(opp)

%tt_1 = [1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 1 1 1 1 0 1];
%tt_2 = [1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1];

%tt_1 = [1 1 1 0 0 1 1 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 1 1 1 0 0 1];
%tt_2 = [1 0 1 0 0 1 1 0 0 1 0 1 1 1 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0];

with_off = 1;

if (with_off)

    went_off_1 = false;
    went_off_2 = false;
    new_frame_2 = false;
    quant_frame_2 = 1;
    contacts_off = 0;
    counted = false;
    contacts=[];
    
    for i=1:1:length(tt_1)
        if(tt_1(i)==1)
            if (went_off_1)
                if (contacts_off~=0)
                    contacts = [contacts contacts_off];
                else
                    contacts = [contacts 0];
                    %contacts = [contacts -quant_frame_2];
                end
                counted = false;
                quant_frame_2 = 1;
                contacts_off = 0;
                if (tt_2(i)==0)
                   went_off_2 = true; 
                else
                   went_off_2 = false;
                end
                went_off_1 = false;
            end
        else
            went_off_1 = true;   
        end
        
        if (tt_2(i)==1)
            if (went_off_2)
                went_off_2 = false;
                quant_frame_2 = quant_frame_2 + 1;
                counted = false;
            end
        else
            went_off_2 = true;
            if (went_off_1 && ~counted)
                contacts_off = contacts_off + 1;
                counted = true;
            end
        end
    end
    
    disp(horzcat('Mean of contacts considering zero contacts = ', num2str(mean(contacts))));
    
    contacts_without_off = contacts(contacts>0);
    disp(horzcat('Mean of contacts without considering zero contacts = ', num2str(mean(contacts_without_off))));
    
    accum=0;
    for i=1:1:max(contacts)
        perc_contact = length(find(contacts==i))/(length(contacts));
        accum = accum + perc_contact;
        disp(horzcat('Ammount (%) of PU1 frames with ', num2str(i),' OFF contacts: ',num2str(perc_contact)));
    end
    
    disp(horzcat('Accum (%) of contact =', num2str(accum)));

else
    
    went_off_1 = false;
    went_off_2 = false;
    new_frame_2 = false;
    quant_frame_2 = 1;
    contacts=[];
    for i=1:1:length(tt_1)
        if(tt_1(i)==1)
            if (went_off_1)
                contacts = [contacts quant_frame_2];
                quant_frame_2 = 1; 
                if (tt_2(i)==0)
                   went_off_2 = true; 
                else
                   went_off_2 = false;
                end
                went_off_1 = false;
            end
        else
            went_off_1 = true;   
        end

        if (tt_2(i)==1)
            if (went_off_2)
                went_off_2 = false;
                quant_frame_2 = quant_frame_2 + 1;
            end
        else
            went_off_2 = true;
        end
    end
    
    disp(horzcat('Mean of contacts  = ', num2str(mean(contacts))));

    accum=0;
    for i=1:1:max(contacts)
        perc_contact = length(find(contacts==i))/length(contacts);
        accum = accum + perc_contact;
        disp(horzcat('Ammount (%) of PU1 frames with ', num2str(i),' PU2 contacts: ',num2str(perc_contact)));
    end
    
    disp(horzcat('Accum (%) of contact =', num2str(accum)));

end
