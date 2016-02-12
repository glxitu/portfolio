%% Simulador de uma Rede Cognitiva descentralizada com "multi channel", a operar utilizando o protocolo MAC proposto - MC-C2RMAC

function [goodput_sum, tx_sim, service_sim, stats_queue] = simulate_multichannel_dualradio_v6(frames, n, ch, idle, cw1, mean_duration_packet, mean_duration_on, lambda)

%% Packet generation. The transmitting SUs and the channel are chosen in order to give priority to SUs with higher backlog in the channels with more availability. Geometric Distributions

% Packet Generation - 0) all to all; 1) half to half; 2) exclusive
all_receive = 2;

%% Probabilities for PU activity
pu_on = (1-idle);                   % PU ON = 1 - PU OFF
%mean_duration_on  = 3;             % Mean of frames on 1 or 7
geo = 1;

pu = [];
if (geo == 1)
    if (idle ~= 1)
        for k=1:1:ch
            mean_duration_off = (1 - pu_on(k)) / pu_on(k)*mean_duration_on;       % Mean of frames off
            l = 0;
            pu_aux = [];
            while (l < frames),                                            % Do the activity map
                on = geornd(1/mean_duration_on)+1;
                ON = ones(1,ceil(on));
                off = geornd(1/mean_duration_off)+1;
                OFF = zeros(1,ceil(off));
                pu_aux = [pu_aux ON OFF];
                l = length(pu_aux);
            end
            pu_aux = pu_aux(1:frames);
            pu = [pu; pu_aux]; 
        end

        %disp(strcat('PU SU ON probabilty EXP = ', num2str(mean(pu(1,:)))));
    else
        pu = zeros(ch, frames);
    end
    
    %disp(strcat('PU ON probabilty EXP = ', num2str(mean(pu(1,:)))));
else
    for i=1:1:ch
        pu_aux = rand(1, frames);
        pu_aux2 = zeros(1, frames);
        pu_aux2(find(pu_aux < pu_on(i))) = 1;
        pu = [pu; pu_aux2];
        %disp(strcat('PU ON probabilty UNI = ', num2str(mean(pu(i,:)))));
    end
end

%% Packet Generation
rarrival = geornd(1/lambda,n/2,2e5) + 1;
packet_arrival=[];
for ii=1:1:n/2
    packet_arrival_aux =[];
    for i=1:length(rarrival(ii,:))
        if rarrival(ii,i) == 1
            packet_arrival_aux = [packet_arrival_aux 1];
        else
            packet_arrival_aux = [packet_arrival_aux zeros(1,rarrival(ii,i)-1) 1];
        end
    end
    if (length(packet_arrival) > length(packet_arrival_aux) && length(packet_arrival) ~= 0)
        packet_arrival = packet_arrival(:,1:length(packet_arrival_aux));
    elseif (length(packet_arrival) < length(packet_arrival_aux) && length(packet_arrival) ~= 0)
        packet_arrival_aux = packet_arrival_aux(1:length(packet_arrival));
    end
    packet_arrival = [packet_arrival; packet_arrival_aux];
end

%% Probabilities of Detection and False Alarm
frame_period = 20e-3;                   % Frame Length 20 ms 
band = 10e3;                            % Bandwidth 10 KHz
sensing_period = 1/(2*band);            % Duration of each sensing sample
number_of_samples = 22;                 % Optimal number of sensing samples
sensing_percentage =  number_of_samples * sensing_period / frame_period;

pd_num = 1;
pfa_num = 0;

pd_aux = rand(ch,frames);
pfa_aux = rand(ch,frames);
pd = ones(ch,frames);
pfa = zeros(ch,frames);
pd(find(pd_aux > pd_num)) = 0;
mean(pd);
pfa(find(pfa_aux < pfa_num)) = 1;
mean(pfa);

%% Channel Availability
for j=1:1:frames 
    slot_idle(:,j) = (1-pd(:,j)).*pu(:,j)+(1-pfa(:,j)).*(1-pu(:,j));
end

%% Initialization
% Protocol variables

channel_state(1:ch) = 0;        % 0 - Idle; 1 - Just assigned; 2 - In use
channel_node_tx(1:ch) = 0;      % ID of the SU transmitting in this channel

state(1:n) = 0;                 % 0 - Idle; 1 - With Packet; 2 - Transmitting; 3 - Receiving 
node(n).packet_sizes = [];      % Number of packets in the queue
node(n).packet_destinies = [];  % ID of the packet's destiny
node_mini_slot(1:n) = 0;        % Selected mini-slot for the RTS
node_transmitted(1:n) = 0;      % Number of packets transmitted
node_backlog(1:n) = 0;

mini_slot(1:cw1) = 0;           % Counts the mini-packets for each RTS phase

access(1:ch,1:frames) = 0;      % Counts the access to the channel

% Stats variables
markov_state(1:n,1:3) = 0;          % 0 - Idle; 1 - Competing; 2 - Transmitting
marked(1:n) = 0;
trans11 = 0;                    % IDLE -> IDLE
trans12 = 0;                    % IDLE -> COMP
trans13 = 0;                    % IDLE -> TRANS NAO EXISTE!
trans21 = 0;                    % COMP -> IDLE
trans22 = 0;                    % COMP -> COMP
trans23 = 0;                    % COMP -> TRANS
trans31 = 0;                    % TRANS -> IDLE NAO EXISTE!
trans32 = 0;                    % TRANS -> COMP
trans33 = 0;                    % TRANS -> TRANS
last = 1;                       % Last state for stats 
node(n).trans_time = [];
node(n).service_time = [];
node(n).init_time = 0;

goodput(1:ch) = 0;
collision(1:ch) = 0;
interference(1:ch) = 0;
unused(1:ch) = 0;
busy(1:ch) = 0;

channels_idle(1:j) = 0;         % Counts idle - without SUs - channels
packet_size = [];
n_sending_rts = [];
psucc = [];
competiu_canal = 0;
passou_canal = [];
competiu_rts = 0;
passou_rts = [];

stats_queue=[];
queue_length=[];
comecou_on = 0;
comecou_off = 0;


%% Protocol
for j=1:1:frames
    
    % Packet Generation
    for i=1:1:n/2
        if packet_arrival(i,j) == 1
            node(i).packet_sizes = [node(i).packet_sizes geornd(1/mean_duration_packet,1,1)+1];
            node(i).packet_destinies = [node(i).packet_destinies ((n/2)+i)];
        end
    end
    
    % Stats
    marked(1:n) = 0;
    channels_idle(j) = length(find(channel_state == 0));
    
    if (length(node(1).packet_sizes) >= 1)
        stats_queue = [stats_queue 1];
        queue_length = [queue_length length(node(1).packet_sizes)];
    else
        stats_queue = [stats_queue 0];
        queue_length = [queue_length 0];
    end
    % end Stats

    % Common Control Channel - it will only work if there is any data channel free
    if (length(find(channel_state == 0)) >= 1)
        % Generate new RTS messages - random destiny, packet length and mini-slot
        node_mini_slot(1:n) = 0;
        mini_slot(1:cw1) = 0;
        for i=1:1:n,
            if (state(i) == 0 && length(node(i).packet_destinies) > 0)          % Only generate RTS for idle SUs and with packet
                node_mini_slot(i) = unidrnd(cw1);
                mini_slot(node_mini_slot(i)) = mini_slot(node_mini_slot(i)) + 1;
                
                state(i) = 1;
            end
        end
        
        % Stats
        if (state(1) == 1)
            competiu_rts = 1;
        end
        n_sending_rts = [n_sending_rts length(find(state==1))];
        psucc = [psucc length(find(mini_slot==1))/length(mini_slot)];
        
        for i=1:1:n
            if (state(i) == 1)
                if (i==1 && last==1)
                    trans12 = trans12 + 1;
                    last = 2;
                elseif (i==1 && last==2)
                    trans22 = trans22 + 1;
                    last = 2;    
                elseif (i==1 && last==3)
                    trans32 = trans32 + 1;
                    last = 2;
                end
                markov_state(i,2) = markov_state(i,2) + 1;
                marked(i) = 1;
            elseif (state(i) == 0)
                if (i==1 && last==1)
                    trans11 = trans11 + 1;
                    last = 1;
                elseif (i==1 && last==2)
                    trans21 = trans21 + 1;
                    last = 1;    
                elseif (i==1 && last==3)
                    trans31 = trans31 + 1;
                    last = 1;
                end
                markov_state(i,1) = markov_state(i,1) + 1;
                marked(i) = 1;
            end
        end
        % end Stats
        
        % Process RTS messages - Ignore two or more RTS packets in the same mini-slot, and busy receivers
        if (length(find(state == 1)) > 0)
            for i=1:1:n,                    
                if (state(i) == 1)              
                    if (mini_slot(node_mini_slot(i)) ~= 1 || state(node(i).packet_destinies(1)) == 2 || state(node(i).packet_destinies(1)) == 3)
                        state(i) = 0;
                    end
                end
            end
        end
        
        % Stats
        if (competiu_rts == 1 && state(1) == 0)
            passou_rts = [passou_rts 0];
            competiu_rts = 0;
        elseif (competiu_rts == 1 && state(1) ==1)
            passou_rts = [passou_rts 1];
            competiu_rts = 0;
        end
        % end Stats
        
        % Update Backlog size
        for i=1:1:n,
            node_backlog(i) = sum(node(i).packet_sizes);
        end
        %disp(strcat('BL:', num2str(node_backlog)));
        
        % Select which SUs will receive a channel, and which channel
        if (length(find(state == 1)) > 0)
            
            % Sort SUs by backlog size - highest to low
            index_sus_competing = find(state==1);
            
            backlog_vec = node_backlog(index_sus_competing);
            backlog_vec = unique(backlog_vec);
            backlog_vec = sort(backlog_vec, 'descend');

            index_sus_sorted = [];
            for m=1:1:length(backlog_vec)
                index_sus = [];
                for mm=1:1:length(index_sus_competing)
                    if (node_backlog(index_sus_competing(mm)) == backlog_vec(m))
                        index_sus = [index_sus index_sus_competing(mm)];
                    end
                end
                index_sus_sorted = [index_sus_sorted index_sus];
            end
            
            % Sort Channels by idle probability - bigger to smaller
            index_channels_free = find(channel_state == 0);
            
            channel_idle_vect = idle(index_channels_free);
            channel_idle_vect = unique(channel_idle_vect);
            channel_idle_vect = sort(channel_idle_vect, 'descend');
            
            index_channels_sorted = [];
            for m=1:1:length(channel_idle_vect)
                index_channel = [];
                for mm=1:1:length(index_channels_free)
                    if ( idle(index_channels_free(mm)) == channel_idle_vect(m))
                        index_channel = [index_channel index_channels_free(mm)];
                    end
                end
                index_channels_sorted = [index_channels_sorted index_channel];
            end
            
            channels_assigned = 0;
            
            % Stats
            if (state(1)==1)
                competiu_canal=1;
            end
            % end Stats

            for i=1:1:length(index_sus_sorted)
                if ((state(index_sus_sorted(i)) ~= 3) && (state(node(index_sus_sorted(i)).packet_destinies(1)) ~= 2) && (state(node(index_sus_sorted(i)).packet_destinies(1)) ~= 3) && (channels_assigned < length(index_channels_sorted)))   % Check if destiny has become busy or if the sender has become a destiny
   
                    state(node(index_sus_sorted(i)).packet_destinies(1)) = 3;       % Rx is set to receive
                    state(index_sus_sorted(i)) = 2;                                 % Tx is set to transmit
                    
                    channels_assigned = channels_assigned + 1;
                    
                    channel_node_tx(index_channels_sorted(channels_assigned)) = index_sus_sorted(i);
                    channel_state(index_channels_sorted(channels_assigned)) = 1;   % Channel is set to "Just Assigned"
                
                    % Stats
                    packet_size = [packet_size node(index_sus_sorted(i)).packet_sizes(1)];
                    node(index_sus_sorted(i)).init_time = j;
                    
                    if (index_sus_sorted(i) == 1 && slot_idle(index_channels_sorted(channels_assigned), j+1) == 1)
                        comecou_off = comecou_off + 1;
                    elseif (index_sus_sorted(i) == 1 && slot_idle(index_channels_sorted(channels_assigned), j+1) == 0)
                        comecou_on = comecou_on + 1;
                    end
                    % end Stats
                    
                elseif (state(index_sus_sorted(i)) ~= 3)
                    
                    state(index_sus_sorted(i)) = 0;
                end
            end
        end
    end
    
    % Stats
    if (competiu_canal==1 && state(1)==2)
        passou_canal = [passou_canal 1];
        competiu_canal = 0;
    elseif (competiu_canal==1)
        passou_canal = [passou_canal 0];
        competiu_canal = 0;
    end
    
    for i=1:1:n
        if (marked(i) == 0)
            if (state(i) == 0)
                if (i==1 && last==1)
                    trans11 = trans11 + 1;
                    last = 1;
                elseif (i==1 && last==2)
                    trans21 = trans21 + 1;
                    last = 1;    
                elseif (i==1 && last==3)
                    trans31 = trans31 + 1;
                    last = 1;
                end
                markov_state(i,1) = markov_state(i,1) + 1;
            elseif (state(i) == 1)
                if (i==1 && last==1)
                    trans12 = trans12 + 1;
                    last = 2;
                elseif (i==1 && last==2)
                    trans22 = trans22 + 1;
                    last = 2;
                elseif (i==1 && last==3)
                    trans32 = trans32 + 1;
                    last = 2;
                end
                markov_state(i,2) = markov_state(i,2) + 1;
            else
                if (i==1 && last==1)
                    trans13 = trans13 + 1;
                    last = 3;
                elseif (i==1 && last==2)
                    trans23 = trans23 + 1;
                    last = 3;
                elseif (i==1 && last==3)
                    trans33 = trans33 + 1;
                    last = 3;
                end
                markov_state(i,3) = markov_state(i,3) + 1;
            end
        end
    end
    % end Stats
    
    % Data Channels
    for k=1:1:ch
        if (channel_state(k) == 2)
            if (slot_idle(k, j) == 1)
                access(k,j) = access(k,j) + 1;
                
                node_transmitted(channel_node_tx(k)) = node_transmitted(channel_node_tx(k)) + 1;
                node(channel_node_tx(k)).packet_sizes(1) = node(channel_node_tx(k)).packet_sizes(1) - 1;
                if (node(channel_node_tx(k)).packet_sizes(1) == 0)
                    state(node(channel_node_tx(k)).packet_destinies(1)) = 0;      % Reset Rx
                    
                    node(channel_node_tx(k)).packet_sizes = node(channel_node_tx(k)).packet_sizes(2:end);
                    node(channel_node_tx(k)).packet_destinies = node(channel_node_tx(k)).packet_destinies(2:end);
                    
                    state(channel_node_tx(k)) = 0;                      % Reset Tx
                    
                    % Stats
                    node(channel_node_tx(k)).service_time = [node(channel_node_tx(k)).service_time j-node(channel_node_tx(k)).init_time];
                    % end Stats
                    
                    channel_state(k) = 0;                               % Free the Channel
                    channel_node_tx(k) = 0;                             % Free the Channel
                end
            end
        elseif (channel_state(k) == 1)
            channel_state(k) = 2;
        end
    end
    
    
    % Stats
    for k=1:1:ch
        
        % Goodput per channel - only one SU and no PUs
        if (access(k,j) == 1 && pu(k,j) == 0)
            goodput(k) = goodput(k) + 1;
        end
        
        % Collision per channel - more than on SU and no PUs
        if (access(k,j) > 1 && pu(k,j) == 0)
            collision(k) = collision(k) + 1;
        end
        
        % Inteference per channel - more than on SU and no PUs
        if (access(k,j) >= 1 && pu(k,j) == 1)
            interference(k) = interference(k) + 1;
        end
        
        % Unused per channel - no SUs and no PUs
        if (access(k,j) == 0 && pu(k,j) == 0)
            unused(k) = unused(k) + 1;
        end
        
        % Busy - Pus
        if (pu(k,j) == 1)
            busy(k) = busy(k) + 1;
        end
    end
end

% Process Stats
if (all_receive == 0)
    for kk=1:1:n
        %for i=2:1:length(node(kk).trans_time)
        %    node(kk).service_time = [node(kk).service_time node(kk).trans_time(i)-node(kk).trans_time(i-1)];
        %end
        node(kk).mean_service_time = mean(node(kk).service_time);
    end
else
    for kk=1:1:(n/2)
        %for i=2:1:length(node(kk).trans_time)
        %    node(kk).service_time = [node(kk).service_time node(kk).trans_time(i)-node(kk).trans_time(i-1)];
        %end
        node(kk).mean_service_time = mean(node(kk).service_time);
    end
end

%goodput_norm = goodput*(1-sensing_percentage)/frames;
goodput_norm = goodput/frames;
collision_norm = collision/frames;
interference_norm = interference/frames;
unused_norm = unused/frames;
busy_norm = busy/frames;

valid_norm = goodput_norm + collision_norm + unused_norm + busy_norm;

node_transmitted;

disp(strcat('Goodput Sim      = ', num2str(goodput_norm,4)));
disp(strcat('Collision Sim    = ', num2str(collision_norm,4)));
disp(strcat('Interference Sim = ', num2str(interference_norm,4)));
disp(strcat('Unused Sim       = ', num2str(unused_norm,4)));
disp(strcat('Busy Sim         = ', num2str(busy_norm,4)));
disp(' ');
for kk=1:1:(n/2)
    disp(strcat('Markov State SU ',num2str(kk),' = ', num2str(markov_state(kk,:)/frames,4)));
end
disp(' ');
for kk=1:1:(n/2)
    disp(strcat('Mean Service Time SU ',num2str(kk),' = ', num2str(node(kk).mean_service_time,4)));
end

pcanal_livre = 1-length(find(channels_idle==0))/length(channels_idle);
psucc_rts = length(find(passou_rts==1))/length(passou_rts);
psucc_canal = length(find(passou_canal==1))/length(passou_canal);

disp(' ');
disp(strcat('Mean sending RTS:', num2str(mean(n_sending_rts))));
disp(strcat('Prob Channel Free:', num2str(pcanal_livre)));
disp(strcat('Prob Success RTS:', num2str(psucc_rts)));
disp(strcat('Prob Channel Ass:', num2str(psucc_canal)));
disp(' ');
disp(strcat('IDLE->IDLE:', num2str(trans11/markov_state(1,1))));
disp(strcat('IDLE->COMP:', num2str(trans12/markov_state(1,1))));
disp(strcat('IDLE->TRANS:', num2str(trans13/markov_state(1,1))));
disp(strcat('COMP->IDLE:', num2str(trans21/markov_state(1,2))));
disp(strcat('COMP->COMP:', num2str(trans22/markov_state(1,2))));
disp(strcat('COMP->TRANS:', num2str(trans23/markov_state(1,2))));
disp(strcat('TRANS->IDLE:', num2str(trans31/markov_state(1,3))));
disp(strcat('TRANS->COMP:', num2str(trans32/markov_state(1,3))));
disp(strcat('TRANS->TRANS:', num2str(trans33/markov_state(1,3))));
disp(' ');
disp(strcat('Prob Packet Queue:', num2str(mean(stats_queue))));
disp(strcat('Mean Packets Queue:', num2str(mean(queue_length))));
disp(strcat('Start Off:', num2str(comecou_off/(comecou_off+comecou_on))));

% Stats to generate figures
goodput_sum = sum(goodput_norm);
tx_sim = mean(markov_state(1:(n/2),3)/frames);
service_sim = 0;
for i=1:(n/2)
    service_sim = service_sim + node(i).mean_service_time;
end
service_sim = service_sim / i;

%% PMF Channels without SUs
% [ff,xx] = hist(channels_idle,0:max(channels_idle));
% figure
% bar(xx,ff/length(channels_idle),'b');
% xlabel('# Idle Channels');ylabel('Probability');
% grid;
%% PMF Packet Size
[ff,xx] = hist(packet_size,0:max(packet_size));
figure
bar(xx,ff/length(packet_size),'b');
xlabel('Packet Size');ylabel('Probability');
grid;
mean(packet_size)
%% Packet Distribution
% [ff,xx] = hist(queue_length,0:max(queue_length));
% figure
% bar(xx,ff/length(queue_length),'b');
% xlabel('# Packets in the Queue');ylabel('Probability');
% mean(queue_length);
% grid;
%% PMF Service Time
[ff,xx] = hist(node(1).service_time,0:max(node(1).service_time));
figure
bar(xx,ff/length(node(1).service_time),'r');
xlabel('Service Time');ylabel('Probability');
grid;
mean(node(1).mean_service_time)
var(node(1).service_time)

