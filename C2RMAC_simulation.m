%% Simulador de uma Rede Cognitiva descentralizada com "single channel", a operar utilizando o protocolo MAC proposto - C2RMAC

function [throughput_sim, goodput_sim_entre_sus, goodput_sim, psucc, pgood, steady_ap, steady_su] = simulate_n_nodes_double_window_no_idle_frames_v2(n, nframes, cw1, cw2, idle_ap, idle_su, lambda)

disp('-------------------------------');
disp(strcat('Number of nodes = ', num2str(n)));

%% Probabilities for PU activity
pu_on_ap = (1-idle_ap);
pu_on_su = (1-idle_su); 
mean_duration_on  = 7;           % Mean of frames on 1 or 7

exp = 0;
same_activity_vector_as_ap = 1;
same_activity_vector_for_sus = 1;

%% Compute PU's activity array - AP
if (exp == 1)
    if (idle_ap ~= 1)
        mean_duration_off = (1 - pu_on_ap) / pu_on_ap*mean_duration_on;       % Mean of frames off
        k = 0;
        pu_ap = [];
        while (k < nframes),                                          % Do the activity map
            on = exprnd(mean_duration_on);
            ON = ones(1,ceil(on));
            off = exprnd(mean_duration_off);
            OFF = zeros(1,ceil(off));
            pu_ap = [pu_ap ON OFF];
            k = length(pu_ap);
        end
    else
        pu_ap = zeros(1, nframes);
    end
    
    disp(strcat('PU AP ON probabilty EXP = ', num2str(mean(pu_ap))));
else
    pu_ap_aux = rand(1, nframes);
    pu_ap = zeros(1, nframes);
    pu_ap(find(pu_ap_aux < pu_on_ap)) = 1;
    
    disp(strcat('PU AP ON probabilty UNI = ', num2str(mean(pu_ap))));
end


%% Compute PU's activity array - SU
pu_su = [];

if (same_activity_vector_as_ap == 0)
    if (same_activity_vector_for_sus == 0)
        if (exp == 1)
            for i=1:1:n
                if (idle_su ~= 1)
                    mean_duration_off = (1 - pu_on_su) / pu_on_su*mean_duration_on;       % Mean of frames off
                    k = 0;
                    pu_su_aux = [];
                    while (k < nframes),                                          % Do the activity map
                        on = exprnd(mean_duration_on);
                        ON = ones(1,ceil(on));
                        off = exprnd(mean_duration_off);
                        OFF = zeros(1,ceil(off));
                        pu_su_aux = [pu_su_aux ON OFF];
                        k = length(pu_su_aux);
                    end
                else
                    pu_su_aux = zeros(1, nframes);
                end
                pu_su_aux = pu_su_aux(1:nframes);
                pu_su = [pu_su; pu_su_aux]; 
            end

            disp(strcat('PU SU ON probabilty EXP = ', num2str(mean(pu_su(1,:)))));
        else
            pu_su_aux = rand(n, nframes);
            pu_su = zeros(n, nframes);
            pu_su(find(pu_su_aux < pu_on_su)) = 1;

            disp(strcat('PU SU ON probabilty UNI = ', num2str(mean(pu_su(1,:)))));
        end
    else
        if (exp == 1)
            if (idle_su ~= 1)
                mean_duration_off = (1 - pu_on_su) / pu_on_su*mean_duration_on;       % Mean of frames off
                k = 0;
                pu_su_aux = [];
                while (k < nframes),                                          % Do the activity map
                    on = exprnd(mean_duration_on);
                    ON = ones(1,ceil(on));
                    off = exprnd(mean_duration_off);
                    OFF = zeros(1,ceil(off));
                    pu_su_aux = [pu_su_aux ON OFF];
                    k = length(pu_su_aux);
                end
            else
                pu_su_aux = zeros(1, nframes);
            end
            pu_su_aux = pu_su_aux(1:nframes);
            
            for i=1:1:n
                pu_su = [pu_su; pu_su_aux]; 
            end

            disp(strcat('PU SU ON probabilty EXP = ', num2str(mean(pu_su(1,:)))));
        else
            pu_su_aux = rand(1, nframes);
            pu_su_aux2 = zeros(1, nframes);
            pu_su_aux2(find(pu_su_aux < pu_on_su)) = 1;
            
            for i=1:1:n
               pu_su = [pu_su; pu_su_aux2]; 
            end

            disp(strcat('PU SU ON probabilty UNI = ', num2str(mean(pu_su(1,:)))));
        end
    end
else
    for i=1:1:n
        pu_su = [pu_su; pu_ap];
    end
end

%% Estatistica de semelhanca de actividade
% for i=1:1:n
%    disp(strcat('Similiarity AP to SU ',num2str(i),' :', num2str(length(find(pu_su(i,:)==pu_ap))/nframes)));
% end

%% Probabilities of detection and false alarm
frame_period = 20e-3;              % 20 ms ou 100 ms   ML
band = 10e3;                        % 10 KHz
sensing_period = 1/(2*band);        % 25 microSeconds - a amostragem ? realizada a 2 vezes a banda, ou seja duas vezes acima da freq de nyquist
number_of_samples = 22;             
sensing_percentage =  number_of_samples * sensing_period / frame_period;

pd_num_ap = 1;      % 0.914
pd_num_ind = 1;     % 0.914

pfa_num_ap = 0;    % 0.0068
pfa_num_ind = 0;   % 0.0068

if (same_activity_vector_as_ap == 0)
    pd_aux_ind = rand(n, nframes);
    pfa_aux_ind = rand(n, nframes);
    pd_ind = ones(n, nframes);              % 1 indica que detecta com succe; 0 indica que n?o detecta, misdetection
    pfa_ind = zeros(n, nframes);            % 0 indica que nao ocorre falso alarme; 1 indica falso alarme
    pd_ind(find(pd_aux_ind > pd_num_ind)) = 0;
    mean(pd_ind);
    pfa_ind(find(pfa_aux_ind < pfa_num_ind)) = 1;
    mean(pfa_ind);

    pd_aux_ap = rand(1,nframes);
    pfa_aux_ap = rand(1,nframes);
    pd_ap = ones(1,nframes);
    pfa_ap = zeros(1,nframes);
    pd_ap(find(pd_aux_ap > pd_num_ap)) = 0;
    mean(pd_ap);
    pfa_ap(find(pfa_aux_ap < pfa_num_ap)) = 1;
    mean(pfa_ap);
else
    pd_aux_ap = rand(1,nframes);
    pfa_aux_ap = rand(1,nframes);
    pd_ap = ones(1,nframes);
    pfa_ap = zeros(1,nframes);
    pd_ap(find(pd_aux_ap > pd_num_ap)) = 0;
    mean(pd_ap);
    pfa_ap(find(pfa_aux_ap < pfa_num_ap)) = 1;
    mean(pfa_ap);
    
    for i=1:1:n
        pd_ind(i,:) = pd_ap;
        pfa_ind(i,:) = pfa_ap;
    end
end

%% Final statistics
throughput_sim = 0;
collision_sim = 0;
interference_sim = 0;
goodput_sim = 0;
goodput_sim_entre_sus = 0;
no_tx_sim = 0;
stages_int_sim = 0;

packet(1:1:n) = 1;

service_init = 0;
service_time = [];
service_counted = 0;

stage_ap = 0;           % 0 - idle; 1 - emite tone 1; 2 - emite tone 2; 3 - idle para TX;
stage(1:1:n) = 0;       % 0 - idle; 1 - verifica passagem a stage 2; 2 - stage 2; 3 - transmitir;

for i=1:1:n
    node(i).mini_slot_frame = [];
end

frames_to_tx(1:1:n) = 0;                % Transmite passadas x tramas idle (PU)
frames_to_t1 = 0;                       % AP volta a T1 passadas x tramas idle (PU)

access(1:1:nframes) = 0;                % Numero de acessos na trama k

access_mini_slot = 0;

list_compete_cw1 = [];
list_compete_cw2 = [];

% STATS
participaram_cw1(1:1:nframes) = -1;
participaram_cw2(1:1:nframes) = -1;

estado_markov_ap(1:1:5) = 0;
estado_markov_su(1:1:6) = 0;

trans11 = 0;
trans12 = 0;
trans21 = 0;
trans23 = 0;
trans24 = 0;
trans31 = 0;
trans33 = 0;
trans34 = 0;
trans45 = 0;
trans46 = 0;
trans51 = 0;
trans52 = 0;
trans55 = 0;
trans56 = 0;
trans61 = 0;
trans62 = 0;

SUmaiorAP = 0;
SUcw2 = 0;
SUtx = 0;
SUntx = 0;

init_ap_t1 = 0;
init_ap_t11 = 0;
init_ap_t2 = 0;
was_t11 = 0;
was_t2 = 0;
was_t3 = 0;
tempos_ap_t1t1 = [];
tempos_ap_t11t1 = [];
tempos_ap_t2t1 = [];

init_su_cw1 = 1;
init_su_cw11 = 1;
init_su_cw2 = 0;
init_su_tx = 1;
was_cw11 = 0;
was_cw2 = 0;
was_tx = 0;
tempos_su_t2tx = [];
tempos_su_idle = [];

vect_su_txtx = [];
vect_ap_t2t2 = [];
vect_cw2_busys = [];

sus_ciclo_ap = [];
sus_ciclo_ap_aux = 0;

goodput_entre_sus_por_tx = [];
goodput_entre_sus_por_tx_aux = 0;

idles_em_tx = [];

stat_compete_cw1 = [];
stat_went_to_cw2 = [];
stat_compete_cw2 = [];

packet(1:1:n) = Inf;

for k=1:1:nframes,
    
    % Compute SU's packet probability
    %for i=1:1:n, 
    %    packet(i) = packet(i) + poissrnd(lambda, 1);
    %end
    
    if (service_counted == 0)
        if (packet(1) >= 1 )
            service_init = k;
            service_counted = 1;
        end
    end
    
    % Matriz de indicacao de trama idle/busy por SU
    for i=1:1:n,            
        slot_idle_ind(i,k) = (1-pd_ind(i,k))*pu_su(i,k)+(1-pfa_ind(i,k))*(1-pu_su(i,k));
    end
    
    % Matriz de indicacao de trama idle/busy para o AP
    slot_idle_ap(k) = (1-pd_ap(1,k))*pu_ap(k)+(1-pfa_ap(1,k))*(1-pu_ap(k));
    
    % AP
    if (stage_ap == 0 && slot_idle_ap(k) == 0)                                                  % IDLE -> IDLE
        stage_ap = 0;
        
        % STATS
        estado_markov_ap(1) = estado_markov_ap(1) + 1;
        
    elseif (stage_ap == 0 && slot_idle_ap(k) == 1)                                              % IDLE -> T1
        stage_ap = 1;
        list_compete_cw1 = [];
        
        % STATS
        estado_markov_ap(2) = estado_markov_ap(2) + 1;
        
        participaram_cw1(k) = 0;
        
        if ((stage(1) == 3 || stage(1) == 2) && frames_to_tx(1) ~= 0)
            SUmaiorAP = SUmaiorAP + 1;
        end
        
        if (was_t3 == 1)
            tempos_ap_t2t1 = [tempos_ap_t2t1 k-init_ap_t2];
            was_t3 = 0;
        end
        
        if (was_t2 == 1)
            tempos_ap_t1t1 = [tempos_ap_t1t1 k-init_ap_t1];
            was_t2 = 0;
        end
        init_ap_t1 = k;
        
        if (was_t11 == 1)
           tempos_ap_t11t1 = [tempos_ap_t11t1 k-init_ap_t11];
           was_t11 = 0;
        end   
        
    %elseif (stage_ap == 1 && (length(list_compete_cw1) > 0 || pu_ap(k-1) == 1) && slot_idle_ap(k) == 1)   % T1 -> T2
    elseif (stage_ap == 1 && (length(list_compete_cw1) > 0) && slot_idle_ap(k) == 1)   % T1 -> T2
        stage_ap = 2;
        
        % STATS
        estado_markov_ap(4) = estado_markov_ap(4) + 1;
        
        participaram_cw2(k) = 0;
        
        init_ap_t2 = k;
        was_t2 = 1;
        
    elseif (stage_ap == 1 && length(list_compete_cw1) == 0 && slot_idle_ap(k) == 1)             % T1 -> T1
        stage_ap = 1;
        
        % STATS
        participaram_cw1(k) = 0;
        
        estado_markov_ap(2) = estado_markov_ap(2) + 1;
        init_ap_t1 = k;             % Dado que tem SU a concorrer em T1 !
                
    elseif (stage_ap == 1 && length(list_compete_cw1) == 0 && slot_idle_ap(k) == 0)             % T1 -> IDLE
        stage_ap = 0;
        
        % STATS
        estado_markov_ap(1) = estado_markov_ap(1) + 1;
        
    %elseif (stage_ap == 1 && (length(list_compete_cw1) > 0 || pu_ap(k-1) == 1) && slot_idle_ap(k) == 0)   % T1 -> T1w
    elseif (stage_ap == 1 && (length(list_compete_cw1) > 0) && slot_idle_ap(k) == 0)   % T1 -> T1w
        stage_ap = 11;
        
        % STATS
        estado_markov_ap(3) = estado_markov_ap(3) + 1;
        
        init_ap_t11 = k;
        was_t11 = 1;
        
    elseif (stage_ap == 11 && slot_idle_ap(k) == 1)                                             % T1w -> T2
        stage_ap = 2;
        
        % STATS
        estado_markov_ap(4) = estado_markov_ap(4) + 1;
        
        participaram_cw2(k) = 0;
        
        init_ap_t2 = k;
        was_t2 = 1;
        
    elseif (stage_ap == 11 && slot_idle_ap(k) == 0)                                             % T1w -> T1w
        stage_ap = 11;
        
        % STATS
        estado_markov_ap(3) = estado_markov_ap(3) + 1;
       
    % elseif (stage_ap == 2 && (length(list_compete_cw2) > 0 || pu_ap(k-1) == 1))
    elseif (stage_ap == 2 && (length(list_compete_cw2) > 0))                      % T2 -> Tx
        list_compete_cw2 = [];
        stage_ap = 3;
        
        if (slot_idle_ap(k) == 1)
            frames_to_t1 = frames_to_t1 - 1;
        end
        
        % STATS
        estado_markov_ap(5) = estado_markov_ap(5) + 1; 
        
        was_t3 = 1;
        
        vect_ap_t2t2 = [vect_ap_t2t2 k-1];
        
        goodput_entre_sus_por_tx = [goodput_entre_sus_por_tx goodput_entre_sus_por_tx_aux];
        goodput_entre_sus_por_tx_aux = 0;
        
        if (slot_idle_ap(k) == 1)
            idles_em_tx = [idles_em_tx k-init_ap_t2];
        end
        
    elseif (stage_ap == 2 && length(list_compete_cw2) == 0 && slot_idle_ap(k) == 0)             % T2 -> IDLE
        stage_ap = 0;
      
        % STATS
        estado_markov_ap(1) = estado_markov_ap(1) + 1;
        
    elseif (stage_ap == 2 && length(list_compete_cw2) == 0 && slot_idle_ap(k) == 1)             % T2 -> T1
        stage_ap = 1;
        list_compete_cw1 = [];
        
        % STATS
        estado_markov_ap(2) = estado_markov_ap(2) + 1;
        
        participaram_cw1(k) = 0;
        
        if (was_t11 == 1)
           tempos_ap_t11t1 = [tempos_ap_t11t1 k-init_ap_t11];
           was_t11 = 0;
        end
        
        if (was_t2 == 1)
            tempos_ap_t1t1 = [tempos_ap_t1t1 k-init_ap_t1];
            was_t2 = 0;
        end
        init_ap_t1 = k;
        
    elseif (stage_ap == 3 && frames_to_t1 == 0 && slot_idle_ap(k) == 0)                         % Tx -> IDLE
        stage_ap = 0;
        
        % STATS
        estado_markov_ap(1) = estado_markov_ap(1) + 1;
        
        sus_ciclo_ap = [sus_ciclo_ap sus_ciclo_ap_aux];
        sus_ciclo_ap_aux = 0;
        
    elseif (stage_ap == 3 && frames_to_t1 == 0 && slot_idle_ap(k) == 1)                         % Tx -> T1
        stage_ap = 1;
        list_compete_cw1 = [];
        
        % STATS
        estado_markov_ap(2) = estado_markov_ap(2) + 1;
        
        participaram_cw1(k) = 0;
        
        if ((stage(1) == 3 || stage(1) == 2) && frames_to_tx(1) ~= 0)
            SUmaiorAP = SUmaiorAP + 1;
        end
        
        if (was_t3 == 1)
            tempos_ap_t2t1 = [tempos_ap_t2t1 k-init_ap_t2];
            was_t3 = 0;
        end
        
        if (was_t2 == 1)
            tempos_ap_t1t1 = [tempos_ap_t1t1 k-init_ap_t1];
            was_t2 = 0;
        end
        init_ap_t1 = k;
        
        if (was_t11 == 1)
           tempos_ap_t11t1 = [tempos_ap_t11t1 k-init_ap_t11];
           was_t11 = 0;
        end
        
        sus_ciclo_ap = [sus_ciclo_ap sus_ciclo_ap_aux];
        sus_ciclo_ap_aux = 0;
        
    elseif (stage_ap == 3 && frames_to_t1 > 0)                                                  % Tx -> Tx
        stage_ap = 3;
        
        if (slot_idle_ap(k) == 1)
            frames_to_t1 = frames_to_t1 - 1;
        end
        
        % STATS
        estado_markov_ap(5) = estado_markov_ap(5) + 1;
        
        if (slot_idle_ap(k) == 1)
            idles_em_tx = [idles_em_tx k-init_ap_t2];
        end
      
    end
      
    % SUs
    for i=1:1:n,
        if (stage(i) == 0)              %% ESTAVA IDLE
            if (slot_idle_ind(i,k) == 1 && packet(i) >= 1 && stage_ap == 1)
                
                stage(i) = 1;                        
                node(i).mini_slot_frame = ceil(cw1*rand(1));
                
                access_mini_slot = access_mini_slot + 1;
                
                list_compete_cw1 = [list_compete_cw1 i];
                
                % STATS
                participaram_cw1(k) = participaram_cw1(k) + 1;
                if (i == 1)
                    trans12 = trans12 + 1;
                    estado_markov_su(2) = estado_markov_su(2) + 1;
                    
                    if (was_tx == 1)
                        tempos_su_idle = [tempos_su_idle k-init_su_tx-1];
                        was_tx = 0;
                        was_cw2 = 0;
                        was_cw11 = 0; 
                    elseif(was_tx == 0 && was_cw2 == 1)
                        tempos_su_idle = [tempos_su_idle k-init_su_cw2-1];
                        was_cw2 = 0;
                        was_cw11 = 0;
                    elseif(was_tx == 0 && was_cw2 == 0 && was_cw11 == 1)
                        tempos_su_idle = [tempos_su_idle k-init_su_cw11];
                        was_cw11 = 0;
                    else
                        tempos_su_idle = [tempos_su_idle k-init_su_cw1-1];
                    end
                    
                    init_su_cw1 = k;
                    
                end
            else
               
                stage(i) = 0;
                
                % STATS
                if (i == 1)
                    trans11 = trans11 + 1;
                    estado_markov_su(1) = estado_markov_su(1) + 1;    
                end
            end
        elseif (stage(i) == 1)          %% ESTAVA T1
            if (slot_idle_ind(i,k) == 1 && ismember(i, list_compete_cw1) == 1 && stage_ap == 2)
                
                list_compete_cw1 = list_compete_cw1(list_compete_cw1 ~= i);
                
                stage(i) = 2;                                           
                node(i).mini_slot_frame = ceil(cw2*rand(1));
                
                access_mini_slot = access_mini_slot + 1;
               
                list_compete_cw2 = [list_compete_cw2 i];
                
                % STATS
                participaram_cw2(k) = participaram_cw2(k) + 1;
                if (i == 1)
                    trans24 = trans24 + 1;
                    estado_markov_su(4) = estado_markov_su(4) + 1;
                    
                    SUcw2 = SUcw2 + 1;
                    
                    init_su_cw2 = k;
                    was_cw2 = 1;
                end
                
            elseif (ismember(i, list_compete_cw1) == 1 && stage_ap == 11)
                
                % Retirar daqui e colocar no else seguinte
                list_compete_cw1 = list_compete_cw1(list_compete_cw1 ~= i);
                
                stage(i) = 11;
                
                % STATS
                if (i == 1)
                    trans23 = trans23 + 1;
                    estado_markov_su(3) = estado_markov_su(3) + 1;
                end
                    
            else
                
                % Retirado do elseif anterior
                list_compete_cw1 = list_compete_cw1(list_compete_cw1 ~= i);
                
                stage(i) = 0;
                
                % STATS
                if (i == 1)
                    trans21 = trans21 + 1;
                    estado_markov_su(1) = estado_markov_su(1) + 1;
                end
                
            end
        elseif (stage(i) == 11)          %% ESTAVA T11
            if (slot_idle_ind(i,k) == 1 && stage_ap == 2)
                            
                stage(i) = 2;                                           
                node(i).mini_slot_frame = ceil(cw2*rand(1));
                
                access_mini_slot = access_mini_slot + 1;
                
                list_compete_cw2 = [list_compete_cw2 i];
                                
                % STATS
                participaram_cw2(k) = participaram_cw2(k) + 1;
                if (i == 1)
                    trans34 = trans34 + 1;
                    estado_markov_su(4) = estado_markov_su(4) + 1;
                    
                    SUcw2 = SUcw2 + 1; 
                    
                    init_su_cw2 = k;
                    was_cw2 = 1;
                end
  
            elseif (stage_ap == 11)
                
                stage(i) = 11;
                
                % STATS
                if (i == 1)
                    trans33 = trans33 + 1;
                    estado_markov_su(3) = estado_markov_su(3) + 1;
                end
            
            else
                
                stage(i) = 0;
                
                % STATS
                if (i == 1)
                    trans31 = trans31 + 1;
                    estado_markov_su(1) = estado_markov_su(1) + 1;
                end 
            end
        elseif (stage(i) == 2)          %% ESTAVA T2
 
            stage(i) = 3;

            if (slot_idle_ind(i,k) == 1)
                frames_to_tx(i) = frames_to_tx(i) - 1;
                if (frames_to_tx(i) == 0)
                    access(k) = access(k) + 1;
                    packet(i) = packet(i) - 1;

                    % STATS
                    sus_ciclo_ap_aux = sus_ciclo_ap_aux + 1;
                    
                    if (i == 1)
                        trans46 = trans46 + 1;
                        estado_markov_su(6) = estado_markov_su(6) + 1;
                        vect_su_txtx = [vect_su_txtx k];

                        SUtx = SUtx + 1;

                        tempos_su_t2tx = [tempos_su_t2tx k-init_su_cw2];
                        was_tx = 1;
                        init_su_tx = k;
                        
                        service_time = [service_time k-service_init];
                        service_counted = 0;
                        if (packet(i) >= 1)
                            service_init = k;
                            service_counted = 1;
                        end
                    end
                else
                    if (i == 1)
                        trans45 = trans45 + 1;
                        estado_markov_su(5) = estado_markov_su(5) + 1;
                    end
                end
            else
                if (i == 1)
                    trans45 = trans45 + 1;
                    estado_markov_su(5) = estado_markov_su(5) + 1;
                end
            end
              
        elseif (stage(i) == 3)          %% ESTAVA T3
            if (slot_idle_ind(i,k) == 1 && stage_ap == 1 && packet(i) >= 1)
            
                stage(i) = 1;
                node(i).mini_slot_frame = ceil(cw1*rand(1));
                
                access_mini_slot = access_mini_slot + 1;
                
                list_compete_cw1 = [list_compete_cw1 i];
                
                % STATS
                participaram_cw1(k) = participaram_cw1(k) + 1;
                if (i == 1 && frames_to_tx(1) ~= 0)
                    trans52 = trans52 + 1;
                    estado_markov_su(2) = estado_markov_su(2) + 1;
                    
                    SUntx = SUntx + 1; 
                    
                    tempos_su_idle = [tempos_su_idle k-init_su_cw2-1];
                    was_cw2 = 0;
                    init_su_cw1 = k;
                    
                elseif (i == 1 && frames_to_tx(1) == 0)
                    trans62 = trans62 + 1;
                    estado_markov_su(2) = estado_markov_su(2) + 1;
                    
                    tempos_su_idle = [tempos_su_idle 0];
                    was_cw2 = 0;
                    was_cw11 = 0;
                    was_tx = 0;
                    init_su_cw1 = k;
                end
                
                frames_to_tx(i) = 0;
            
            elseif (slot_idle_ind(i,k) == 0 && stage_ap == 1)
                
                stage(i) = 0;
                
                % STATS
                if (i == 1 && frames_to_tx(1) ~= 0)
                    trans51 = trans51 + 1;
                    estado_markov_su(1) = estado_markov_su(1) + 1;
                    
                    SUntx = SUntx + 1;

                elseif (i == 1 && frames_to_tx(1) == 0)
                    trans61 = trans61 + 1;
                    estado_markov_su(1) = estado_markov_su(1) + 1;
                end
                
                frames_to_tx(i) = 0;
                
            elseif (frames_to_tx(i) == 0 && stage_ap ~= 1)
                
                stage(i) = 0;
                
                % STATS
                if (i == 1)
                    trans61 = trans61 + 1;
                    estado_markov_su(1) = estado_markov_su(1) + 1;
                end
   
            else
                
                stage(i) = 3;
                
                if (slot_idle_ind(i,k) == 1)
                    frames_to_tx(i) = frames_to_tx(i) - 1;
                    if (frames_to_tx(i) == 0)
                        access(k) = access(k) + 1;
                        packet(i) = packet(i) - 1;
                        
                        % STATS
                        sus_ciclo_ap_aux = sus_ciclo_ap_aux + 1;
                        
                        if (i == 1)
                            trans56 = trans56 + 1;
                            estado_markov_su(6) = estado_markov_su(6) + 1;
                            vect_su_txtx = [vect_su_txtx k];
                            
                            SUtx = SUtx + 1;
                            
                            tempos_su_t2tx = [tempos_su_t2tx k-init_su_cw2];                            
                            was_tx = 1;
                            init_su_tx = k;
                            
                            
                           	service_time = [service_time k-service_init];
                            service_counted = 0;
                            if (packet(i) >= 1)
                                service_init = k;
                            service_counted = 1;
                            end
                        end
                    else
                        if (i == 1)
                            trans55 = trans55 + 1;
                            estado_markov_su(5) = estado_markov_su(5) + 1;
                        end
                    end
                else
                    if (i == 1)
                        trans55 = trans55 + 1;
                        estado_markov_su(5) = estado_markov_su(5) + 1;
                    end
                end  
                
            end
        end
    end
     
    if (stage_ap == 1)
        stat_compete_cw1 = [stat_compete_cw1 length(list_compete_cw1)]; 
    end
    if (stage_ap == 2)
        stat_compete_cw2 = [stat_compete_cw2 length(list_compete_cw2)];
    end
         
    if ((stage_ap == 1 || stage_ap == 2) && access_mini_slot ~= 0 )

        cws = [];
        for i=1:1:n,
            cws = [cws node(i).mini_slot_frame];
        end
        cws_u = unique(cws);

        % STATS
        if (stage_ap == 2)
            frames_to_t1 = length(cws_u);
            vect_cw2_busys = [vect_cw2_busys frames_to_t1];
        end

        if (pu_ap(k) == 1)                      % Tem PU
            min_mini_slot = 1;
        else                                    % Nao tem PU
            min_mini_slot = min(cws_u);
        end
        
        for i=1:1:n,
            if (slot_idle_ind(i,k) == 1)
                if (stage(i) == 1)              % Stage 1 - verifica se vai passar ao stage 2 ou ao stage 0

                    if (min(node(i).mini_slot_frame) ~= min_mini_slot)     
                        list_compete_cw1 = list_compete_cw1(list_compete_cw1 ~= i);
                    end
                    node(i).mini_slot_frame = [];
                    
                elseif (stage(i) == 2)          % Stage 2 - calcula quantas tramas IDLE tem de ver antes de transmitir
                    if (min_mini_slot == min(node(i).mini_slot_frame))
                        frames_to_tx(i) = 1;
                    else
                        frames_to_tx(i) = length(find(cws_u < min(node(i).mini_slot_frame))) + 1;
                    end
                    node(i).mini_slot_frame = [];                    
                end   
            end
        end
       
    end
   
    % STATS
    if (stage_ap == 1)
        stat_went_to_cw2 = [stat_went_to_cw2 length(list_compete_cw1)];
    end
    
    % Accesso de 1 ou mais SUs, independentemente do estado do PU no AP      
    if (access(k) >= 1)
        throughput_sim = throughput_sim + 1;
    end
       
    % 1 ou mais SUs colidirem, independentemente do estado do PU no AP
    if (access(k) > 1)                                                  % if ((access(k) > 1 && pu_ap(k) == 0) || (access(k) == 1 && access_mini_slot >= 1) )
        collision_sim = collision_sim + 1;
    end
    
    % 1 ou mais SUs transmitirem e o AP estar ocupado
    if (access(k) >= 1 && pu_ap(k) == 1)                                % if ((access(k) >= 1 && pu_ap(k) == 1) || (access_mini_slot >=1 && pu_ap(k) == 1))    
        interference_sim = interference_sim + 1;
    end  
    
    % 1 ou mais SUs transmitirem, e o AP estar livre (possivel utilizacao de MPR)
    if (access(k) >= 1 && pu_ap(k) == 0)
        goodput_sim_entre_sus = goodput_sim_entre_sus + 1;
        
        goodput_entre_sus_por_tx_aux = goodput_entre_sus_por_tx_aux + 1;
    end
    
    % 1 SU transmitir e o AP estar livre
    if (access(k) == 1 && pu_ap(k) == 0)                                % if (access(k) == 1 && pu_ap(k) == 0 && access_mini_slot == 0)
        goodput_sim = goodput_sim + 1;
    end
    
    % Nenhum SU transmitir e o meio estar livre
    if (access(k) == 0 && pu_ap(k) == 0)
        no_tx_sim = no_tx_sim + 1;
    end
    
    % Trama usada para Stage 1 ou Stage 2
    if (stage_ap == 1 || stage_ap == 2)
        stages_int_sim = stages_int_sim + 1;
    end
    
    access_mini_slot = 0;
    
end

SUcw2;
SUtx;
SUntx;
SUmaiorAP;
ProbAPmaiorSUv1 = SUtx/SUcw2;
ProbAPmaiorSUv2 = 1-(SUmaiorAP/SUcw2); 

disp(strcat('IDLE SU = ', num2str(mean(slot_idle_ind(1,:)))));
disp(strcat('IDLE AP = ', num2str(mean(slot_idle_ap(:)))));

disp(strcat('Pbeta = ', num2str(ProbAPmaiorSUv1)));
disp(strcat('PSU1 = ', num2str(1-length(find(participaram_cw1==0))/estado_markov_ap(2))));
disp(strcat('PSU2 = ', num2str(1-length(find(participaram_cw2==0))/estado_markov_ap(4))));

participaram_cw1(participaram_cw1==-1) = [];
disp(strcat('Participaram CW1 = ', num2str(mean(participaram_cw1))));

participaram_cw2(participaram_cw2==-1) = [];
disp(strcat('Participaram CW2 = ', num2str(mean(participaram_cw2))));

%% Probabilidades de seady state do AP e SU e transicao do SU
steady_ap = estado_markov_ap/nframes
steady_su = estado_markov_su/nframes

% disp(strcat('Trans IDLE->IDLE = ', num2str(trans11/estado_markov_su(1))));
% disp(strcat('Trans IDLE->CW1 = ', num2str(trans12/estado_markov_su(1))));
% 
% disp(strcat('Trans CW1->CW11 = ', num2str(trans23/estado_markov_su(2))));
% disp(strcat('Trans CW1->CW2 = ', num2str(trans24/estado_markov_su(2))));
% disp(strcat('Trans CW1->IDLE = ', num2str(trans21/estado_markov_su(2))));
% 
% disp(strcat('Trans CW11->CW2 = ', num2str(trans34/estado_markov_su(3))));
% disp(strcat('Trans CW11->CW11 = ', num2str(trans33/estado_markov_su(3))));
% disp(strcat('Trans CW11->IDLE = ', num2str(trans31/estado_markov_su(3))));
% 
% disp(strcat('Trans CW2->CW2TOT = ', num2str(trans45/estado_markov_su(4))));
% disp(strcat('Trans CW2->TX = ', num2str(trans46/estado_markov_su(4))));
% 
% disp(strcat('Trans CW2TOT->IDLE = ', num2str(trans51/estado_markov_su(5))));
% disp(strcat('Trans CW2TOT->CW1 = ', num2str(trans52/estado_markov_su(5))));
% disp(strcat('Trans CW2TOT->CW2TOT = ', num2str(trans55/estado_markov_su(5))));
% disp(strcat('Trans CW2TOT->TX = ', num2str(trans56/estado_markov_su(5))));
% 
% disp(strcat('Trans TX->IDLE = ', num2str(trans61/estado_markov_su(6))));
% disp(strcat('Trans TX->CW1 = ', num2str(trans62/estado_markov_su(6))));

%% Gr?ficos de PMFs relacionadas com n1, n2, n2a, cw2b, AP T2->T1, AP T1->T1, SU T2->Tx, SU em IDLE, SU Tx->Tx
% n1
% [fap,xap] = hist(stat_compete_cw1,0:max(stat_compete_cw1));
% 
% mean(stat_compete_cw1)
% 
% figure
% bar(xap,fap/length(stat_compete_cw1),'b');
% grid;

% % n2
[fap,xap] = hist(stat_went_to_cw2,0:max(stat_went_to_cw2));

mean(stat_went_to_cw2)

figure
bar(xap,fap/length(stat_went_to_cw2),'b');
grid;
hold on

n2(1:n)=0;
for k=1:1:n-1
    for i=0:1:cw1-2
        n2(k)=n2(k)+((cw1-i)^(n)/(cw1^(n)))*nchoosek(n, k)*(1/(cw1-i))^(k)*(1-1/(cw1-i))^(n-k);
    end
end
n2(n)=cw1/(cw1)^n;

plot(1:1:n,n2,'r*')

n2theo_mean=0;
for k=1:n
    n2theo_mean=n2theo_mean + k*n2(k);
end

n2theo_mean


% % n2a
% [fap,xap] = hist(participaram_cw2,0:max(participaram_cw2));
% 
% mean(participaram_cw2)
% 
% figure
% bar(xap,fap/length(participaram_cw2),'b'); 
% grid;

% % cw2b
% [fap,xap] = hist(vect_cw2_busys,0:max(vect_cw2_busys));
% 
% mean(vect_cw2_busys)
% 
% figure
% bar(xap,fap/length(vect_cw2_busys),'b');
% grid;

% AP T2 -> T1
% [fap,xap] = hist(tempos_ap_t2t1,1:max(tempos_ap_t2t1));
% xap = [0 xap];
% fap = [0 fap];
% 
% mean(tempos_ap_t2t1)
% 
% figure
% bar(xap,fap/length(tempos_ap_t2t1),'b'); 
% grid;

% % AP T1 -> T1
% [fap,xap] = hist(tempos_ap_t1t1,1:max(tempos_ap_t1t1));
% xap = [0 xap];
% fap = [0 fap];
% 
% mean(tempos_ap_t1t1)
% 
% figure
% bar(xap,fap/length(tempos_ap_t1t1),'b'); 
% grid;

% % SU T2 -> Tx
% [fsu,xsu] = hist(tempos_su_t2tx,1:max(tempos_su_t2tx));
% xsu = [0 xsu];
% fsu = [0 fsu];
% 
% mean(tempos_su_t2tx)
% 
% figure
% bar(xsu,fsu/length(tempos_su_t2tx),'b'); 
% grid;

% SU em IDLE
% [fsu,xsu] = hist(tempos_su_idle,0:max(tempos_su_idle));
% 
% mean(tempos_su_idle)
%  
% figure
% bar(xsu,fsu/length(tempos_su_idle),'b');
% grid;

% SU Tx -> Tx
% [fsu,xsu] = hist(service,1:max(service));
% 
% mean(service);
% 
% [phat,pci] = gamfit(service);
% phat(1)*phat(2);
% 
% figure
% grid on;
% bar(xsu,fsu/length(service),'b'); hold on
% plot(xsu,gampdf(xsu,phat(1),phat(2)),'r*-');
% grid;

%% Outras Estatisticas

idle_su = mean(slot_idle_ind(1,:));
idle_ap = mean(slot_idle_ap(:));

t2_ap = [0];
for i=2:1:length(vect_ap_t2t2)
    t2_ap(i-1) = vect_ap_t2t2(i)-vect_ap_t2t2(i-1);
end
mean(t2_ap);

tx = [0];
service = [0];

for i=2:1:length(vect_su_txtx)
    tx(i-1) = 1/(vect_su_txtx(i)-vect_su_txtx(i-1));
    service(i-1) = (vect_su_txtx(i)-vect_su_txtx(i-1));
end
tx = tx*idle_su;
tx_mean = mean(tx);
tx_std = std(tx);
int_conf_tx = 2.576*(tx_std/sqrt(length(tx)));

service_mean = mean(service_time);
service_std = std(service_time);
service_var = var(service_time);
int_conf_service = 2.576*(service_std/sqrt(length(tx)));

sus_ciclo_ap_mean = mean(sus_ciclo_ap);

goodput_entre_sus_por_tx_mean = mean(goodput_entre_sus_por_tx);
goodput_entre_sus_por_tx_mean/mean(t2_ap);

%% Estatisticas de Throughput e Goodput
throughput_sim = throughput_sim/nframes;
goodput_sim = goodput_sim/nframes;
goodput_sim = goodput_sim*(1-sensing_percentage);
goodput_sim_entre_sus = goodput_sim_entre_sus/(nframes);
collision_sim = collision_sim/nframes;
no_tx_sim = no_tx_sim/nframes;
interference_sim = interference_sim/(nframes);
idle_slots = no_tx_sim/(mean(slot_idle_ap(:)));

pidle = (1-steady_su(6))^n;
psucc = n*steady_su(6)*(1-steady_su(6))^(n-1);
pcol = 1 - pidle - psucc;
pgood = psucc*(mean(slot_idle_ap(:)));
pint = psucc - pgood;

Psu2 = 1-length(find(participaram_cw2==0))/estado_markov_ap(4);
Tciclo_AP = 1/(steady_ap(4)*Psu2);
Tciclo_SU = 1/steady_su(6);

avg_sus_ciclo_ap = (Tciclo_AP/Tciclo_SU)*n;
avg_t2t1 = Tciclo_AP*steady_ap(5);
tau_su = 1/avg_t2t1;
idle_por_ciclo_ap = avg_t2t1*(1-tau_su)^avg_sus_ciclo_ap;
S = (avg_t2t1-idle_por_ciclo_ap)/Tciclo_AP;
G = ((avg_t2t1-idle_por_ciclo_ap)*idle_ap)/Tciclo_AP;

disp(strcat('Throu Sim = ', num2str(throughput_sim)));
disp(strcat('Coll Sim = ', num2str(collision_sim)));
disp(strcat('Int Sim = ', num2str(interference_sim)));
disp(strcat('Idle Sim = ', num2str(idle_slots)));
disp(strcat('Good Good Sim = ', num2str(goodput_sim)));
disp(strcat('Good Sim = ', num2str(goodput_sim_entre_sus)));
disp(' ');
disp(strcat('Service Mean Sim = ', num2str(service_mean)));
disp(strcat('Service Var Sim = ', num2str(service_var)));
disp(strcat('Prob of Queue Empty =', num2str(1-(steady_su(2)/steady_ap(2)))));
disp(' ');
disp(strcat('Throu Theo Prof = ', num2str(psucc)));
disp(strcat('Coll Theo Prof = ', num2str(pcol)));
disp(strcat('Int Theo Prof = ', num2str(pint)));
disp(strcat('Idle Theo Prof = ', num2str(pidle)));
disp(strcat('Good Theo Prof = ', num2str(pgood)));
disp(' ');
disp(strcat('Throu Theo Mig = ', num2str(S)));
disp(strcat('Good Theo Mig = ', num2str(G)));

%[S, G] = plot_distribAPeSU(n, cw1, cw2, idle_ap, idle_su, steady_su, steady_ap, same_activity_vector, tempos_ap_t1t1, tempos_ap_t11t1, tempos_ap_t2t1, idles_em_tx, tempos_su_t2tx, tempos_su_idle, service);
