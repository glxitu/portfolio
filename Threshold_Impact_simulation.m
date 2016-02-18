%% Simulador de uma Rede Cognitiva descentralizada com "single channel" para avalia??o do impacto do treshold de energia a adoptar pelos utilizadores secund?rios

clc;
clear;
close all;

% --- Vari?veis --- %
PU = 70;
ratio = 15;
simu = 1;
% ------------------ %

% --- Amostragem --- %
W = 10*10^3; 
T = 1/(2*W);
sample_time = T; 
% ------------------ %

% --- Dura??o da Frame --- %
T_sense = 21.3 * 10^(-3) * 0.02;
T_data =  21.3 * 10^(-3) * (1 - 0.02);
T_frame = T_sense +  T_data;
disp('T_frame');
disp(T_frame);
% ------------------------ %

for simu=1:1:25
    % --- Comportamento PU --- %
    ficheiro = strcat('comportamento_PU/p',num2str(PU),'_',num2str(ratio),'tf_v',num2str(simu),'.txt');
    aux = load(ficheiro,'-ascii');
    transicoes_meio = [0 cumsum(aux')];
    % ------------------------ %

    % --- Tempo de simula??o --- %
    T_simulation = floor(sum(aux)/T_frame)*T_frame;
    disp('T_simulation');
    disp(T_simulation);
    % -------------------------- %

    % --- Sinal do PU --- %
    SNR = 5;
    signalPU_var = 0;  
    lambda = 10^(SNR/10);
    % ------------------- %

    if(mod(length(transicoes_meio),2)==0)
        time_on = sum(transicoes_meio(3:2:end) - transicoes_meio(2:2:end-1));
        time_off = sum(transicoes_meio(2:2:end) - transicoes_meio(1:2:end));
    else
        time_on = sum(transicoes_meio(3:2:end) - transicoes_meio(2:2:end));
        time_off = sum(transicoes_meio(2:2:end) - transicoes_meio(1:2:end-1));
    end;

    % --- Calculo de p do PU --- %
    p = time_on/(time_on + time_off);
    disp('P_ON:');
    disp(p);
    disp('Tempo real de simula??o : Time_on + Time_off');
    disp(time_on + time_off);
    % -------------------------- %

    % --- Abrir ficheiro de resultados e passo da % de sensing --- %
    passo = 0.0025;
    file_name = strcat('results_simu/res_p',num2str(PU),'_',num2str(ratio),'tf_v',num2str(simu),'.txt');
    results = fopen(file_name,'w');
    Tx = 21.3*10^(-3);
    pos = 0;
    % ------------------------------------------------------------ %

    % Leitura da Interfer?ncia da FU5
    interfFU5 = 0.5;

    for n=0.01:passo:1

        pos = pos + 1;
        % Tempo em segundos
        T_sense = Tx*(n);
        disp('T_sense');
        disp(T_sense);
        T_data = Tx*(1 - n);
        disp('T_data');
        disp(T_data);

        S = floor(T_sense/sample_time);
        disp('S');
        disp(S);

        % -------------------------------------- %
        % ---          Maximiza??es          --- %
        % -------------------------------------- %

        % Maximiza??o pela probabilidade
        %[limiar_Pma, limiar_Pmasucc, limiar_equalpPdPfa, limiar_equalPdPfa, limiar_Pdfa, limiar_Pmad, limiar_TNP]= threshold_computation( SNR, lambda, S, p, interfFU5);
        %threshold = limiar_equalpPdPfa;   %% Limiar de pot?ncia
        filename_thres = strcat('../Fase4/old_param/c4_p',num2str(PU),'.txt');
        threshold_vec = load(filename_thres,'-ascii');
        thres_pos = find(threshold_vec(:,1)==S);
        threshold = threshold_vec(thres_pos,2);
        fprintf(results,'%f\t', threshold);
        disp('Threshold');
        disp(threshold);
        % ---------------------------------------%

        Interf_PU = 0;
        sum_TX = 0;
        num_frameSU = floor(T_simulation/T_frame);
        disp('Numero de Tramas:');
        disp(num_frameSU);
        SU_Detection = zeros([1 floor(num_frameSU)]);
        tramas_utilPu = ones(size(transicoes_meio));

        for j=1:1:num_frameSU,

            %--- slots of sense --- %
                Slot_0   = (j-1) * T_frame + sample_time;
                Slot_end = (j-1) * T_frame + T_sense;
            %---------------------- %

            busy = 0;
            idle = 0;
            % Percorre cada slot de bloco de sensing
            for i=1:1:T_sense/sample_time
                X = Slot_0 + (i-1)*sample_time;
                % Procura ?ltima altera??o do estado do prim?rio
                H = max(find(transicoes_meio <= X));
                % Verifica o estado do prim?rio
                if(mod(H,2)==0)% Existe prim?rio no meio
                    busy =  busy + 1;
                else % N?o exite prim?rio no meio 
                    idle =  idle + 1;
                end;
            end;
            % Faz o sensing do "bloco sentido"
            [result] = make_sensing( idle , busy, SNR, threshold, signalPU_var ,S);
            SU_Detection(j) = result;


            % Para cada trama verifica se o secund?rio colidiu com o prim?rio
            % --- D?bito ?til(colis?o com PU) --- %
            if(SU_Detection(j)==0)

                % --- slots of trama --- %
                % Primeiro slot da trama
                Slot_0 = (j-1) * T_frame;
                % ?ltimo slot da trama 
                Slot_end = j*T_frame;    

                % Para um dado bloco TX verifica se o prim?rio apareceu no meio
                X = max(find(transicoes_meio<= Slot_end));          
                [r,c,v] = find(transicoes_meio > Slot_0+ T_sense & transicoes_meio < Slot_end);

                if(length(c)>1)
                    transicoes_meio(r,c)
                    error('Erro : Mais que uma transi??o numa mesma trama!!!');

                end;

                if(transicoes_meio(X) <= Slot_0 + T_sense)
                    if(mod(X,2)==0) %Tem PU
                        tramas_utilPu(X) = 0;
                        Interf_PU = Interf_PU + (T_frame - T_sense);
                    else % N?o tem PU
                        sum_TX = sum_TX + T_data;
                    end;
                else  
                    if(mod(X,2)==0) % Transi??o de OFF PARA ON a meio da Trama TX
                        tramas_utilPu(X) = 0;
                        K = Slot_end - transicoes_meio(X);
                        if(K<0)
                            K = 0;
                        end;
                        Interf_PU = Interf_PU + K; 

                    else % Transi??o de ON PARA OFF a meio da Trama TX
                        tramas_utilPu(X-1) = 0;
                        K = transicoes_meio(X) - Slot_0;
                        if(K<0)
                            error('Erro : Transi??o de ON PARA OFF a meio da Trama TX');
                        end;
                        Interf_PU = Interf_PU + K;              
                    end;
                end;
            end;
        end;


        % --- D?bito do Prim?rio --- %
        PU_T_slot(pos) = (time_on/T_simulation)*100;
        fprintf(results,'%f \t',PU_T_slot(pos));
        % -------------------------- %

        % --- Goodput PU ao slot --- %
        PU_TU_slot(pos) = ((time_on - Interf_PU)/T_simulation)*100;
        fprintf(results,'%f \t',PU_TU_slot(pos));
        % -------------------------- %

        % --- Goodput PU ? trama --- %
        if(mod(length(transicoes_meio),2)==0)
            time_on_interf = sum((transicoes_meio(3:2:end) - transicoes_meio(2:2:end-1)).*tramas_utilPu(2:2:end-1));
        else
            time_on_interf = sum((transicoes_meio(3:2:end) - transicoes_meio(2:2:end)).*tramas_utilPu(2:2:end));
        end;
        PU_TU_trama(pos) =  (time_on_interf/T_simulation)*100;
        fprintf(results,'%f \t',PU_TU_trama(pos));
        % -------------------------- %

        % --- D?bito do Secund?rio --- %
        SU_T(pos) = sum (SU_Detection) * T_data/T_simulation *100;
        fprintf(results,'%f \t',SU_T(pos));
        % ---------------------------- %

        % --- Goodput SU ao slot --- %
        SU_TU(pos) = sum_TX/T_simulation * 100;
        fprintf(results,'%f \t',SU_TU(pos));
        % -------------------------- %

        % --- Interfer?ncia --- %
        Interferencia_PU(pos) = Interf_PU/T_simulation * 100;
        fprintf(results,'%f \t',Interferencia_PU(pos));
        % --------------------- %

        % --- Percentagem de Sensing --- %
        T_s(pos) = T_sense*100/T_frame;
        fprintf(results,'%f\t',T_s(pos));
        % ------------------------------ %

        % --- Amostras de Sensing --- %
        num_amostras(pos) = S;
        fprintf(results,'%f\t',num_amostras(pos));
        % --------------------------- %

        % --- Tempo de Sensing --- %
        T_sensing(pos) = T_sense;
        fprintf(results,'%f\t',T_sensing(pos));
        % ------------------------ %

        % --- Tempo de Dados --- %
        T_dados(pos) = T_data;
        fprintf(results,'%f\t',T_dados(pos));
        % ---------------------- %

        % --- SNR --- %
        Sim_SNR(pos) = SNR;
        fprintf(results,'%f \t',Sim_SNR(pos));
        % ----------- %

        % --- Percentagem PU Real --- %
        p_pu(pos) = p;
        fprintf(results,'%f \n',p_pu(pos));
        % --------------------------- %

     end;

    % Fecho do ficheiro onde se est?o a guardar os resultados
    fclose(results);

end;
    

