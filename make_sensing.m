function [ result ] = make_sensing(idle , busy, SNR, threshold, signalPU_var ,S)


% --- Sinal do Prim?rio --- %
snr = 10^(SNR/20);
signalPU_m = snr;  %% m?dia da distribui??o normal do sinal
% ------------------------- %

% --- Ru?do --- %
noise_m = 0;     %% m?dia do ru?do (distribui??o Normal) N(0,1)
noise_var = 1;   %% m?dia do ru?do (distribui??o Normal) N(0,1)
% ------------- %

% --- C?lculo da energia --- %
% Calculo do sinal do meio (ru?do + sinal)
% Feita a raiz pois o signalPU_m ? referente ? potencia do sinal

% Ru?do
noise = normrnd(noise_m, sqrt(noise_var),[idle + busy 1]);
% Sinal do PU
signalPU =  normrnd(signalPU_m, sqrt(signalPU_var),[idle + busy 1]);
signalPU(1:idle,1)  = 0;
% Sinal + Ru?do
signal = noise + signalPU;

% Pot?ncia do Sinal do meio
power = sum(((signal.^2)')); %%Coluna com cada posi??o a corresponder ? pot?ncia sentida numa dada trama de sensing
% -------------------------- %

if(power>=threshold)
    result = 1;
else
    result = 0;
end;

end

