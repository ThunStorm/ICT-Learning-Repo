clc;
clear all;
close all;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  The Simulation of an 8-ary Digital Communiation System has been                %
%  implemented with both functions from Communication Toolbox(P3, P4) and  %
%  raw code (P1, P2).                                                                                            %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  


% Define the number of SYMBOLS (equals to 3 bit)
num_symbol = 100000;
% 8-ary QAM
M = 8;
d = 1;
% Gray coded rules (Decimal)
graycode8 = [0 1 3 2 4 5 7 6];
% Define range of EsN0 (symbol/noise)
EsN0 = 0 :  0.1 : 20;
% EsN0 = 0 : 20; % Larger step due to laptop performance limit
% dB value to linear
snr = 10 .^ (EsN0 / 10);
% Random a symbol sequence as original signal (1*num_symbol vector)
msg = randi([0,M-1],1,num_symbol);

%% p3
% Map symbols following gray code and Do modulation
msg_gray = graycode8(msg+1);
msgmod = qammod(msg_gray,M);
% % Constellation diagram
% scatterplot(msgmod);

% %% p1
% % Since functions from Communication toolbox are not suggested to be used
% % A tiring detailed procedure shows below to help us to get a better understanding
% vectors = [-d/2+3*d/2*i, -d/2+d/2*i, -d/2-d/2*i, -d/2-3*d/2*i, d/2-3*d/2*i, d/2-d/2*i, d/2+d/2*i, d/2+3*d/2*i];
% msgmod = zeros(1, num_symbol);
% for i = 1 : num_symbol
%     index_gray = graycode8(msg(i)+1);
%     msgmod(i) = vectors(index_gray+1);
% end

%% p4
% Average energy of each symbol
spow = norm(msgmod) .^ 2 / num_symbol;

% Calculate BER & SER using functions from Communication toolbox
for i = 1 : length(EsN0)
%     Calculate energy of noise
    sigma = sqrt(spow / (2*snr(i)));
%     Blend AWGN in
    r = msgmod + sigma * (randn(1, length(msgmod))+1i * randn(1, length(msgmod)));
%     8-QAM demodulation
    y = qamdemod(r,M);
%     Inverse map messages (demicals obtained)
   decmsg = graycode8(y + 1);
%    Compare signal to the original and calculate BER and SER (A symbol consists of 3 bits)
   [err1, ber(i)] = biterr(msg,decmsg,log2(M));
   [err2, ser(i)] = symerr(msg,decmsg);
end

% % Constellation diagram of the received signal
% scatterplot(r);

% %% p2
% % Since functions from Communication toolbox are not suggested to be used
% % A tiring detailed procedure shows below to help us to get a better understanding
% 
% % Initialization
% Euclidean_distance = zeros(1,M);
% error_count = zeros(1, length(EsN0));
% ser = zeros(1, length(EsN0));
% ber = zeros(1, length(EsN0));
% 
% % Average power of each symbol
% spow = norm(msgmod) .^ 2 / num_symbol;
% 
% for i = 1:length(EsN0)
% %         Calculate energy of noise
%     sigma = sqrt(spow / (2*snr(i)));
% %         Blend AWGN in
%     r = msgmod + sigma * (randn(1, length(msgmod))+1i * randn(1, length(msgmod)));
%     
%     for j = 1 : num_symbol
% %         Calculate Euclidean distances between received symbols and every vector
%         for n = 1 : M
%             Euclidean_distance(n) = norm(r(j) - vectors(n));
%         end
% %         Search the shortest distance, the signal is recognised as the nearest vector
%         pos = find(Euclidean_distance == min(Euclidean_distance));
% %         Count the number of errors
%         if (vectors(pos) ~= msgmod(j))
%             error_count(i) = error_count(i) + 1;
%         end
%     end
% %     SER and BER
%     ser(i) = error_count(i)./num_symbol;
%     ber(i) = error_count(i)./(num_symbol .* log2(M));
% end

%%
% Thereotical SER and BER calculated as follow
ser_theory = 5 / 2 * qfunc(sqrt(snr ./ 3)) - 3 / 2 *(qfunc(sqrt(snr ./ 3))).^2;
ber_theory = 1 / log2(M) * ser_theory;

% Plotting
figure()
subplot(2,1,1)
semilogy(EsN0/log2(M),ber,"o", EsN0/log2(M), ber_theory, "-");
title("8-QAM Bit Error Rate")
xlabel("Eb/N0(dB)");
ylabel("BER");
legend("Result", "Theoretical");
grid on
subplot(2,1,2)
semilogy( EsN0, ser, "*",EsN0, ser_theory, "-");
title("8-QAM Symbol Error Rate")
xlabel("Es/N0(dB)");
ylabel("SER");
legend("Result","Theoretical");
grid on