clc;
clear all;
close all;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  Simulation of a Hamming-coded BPSK System                                               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  

%% Define variables and signal 
n = 7;
k = 4;
% Define the length of data bit
data_length = 400000;
encoded_data_length = data_length / k * n;
% Value assignment to input signal (to generate a bit sequence)
data = round(rand(1, data_length));
% Define coded signal
data_encoded = zeros(1, encoded_data_length);

% Energy per bit
Eb = n / k;
% Eb/N0 on a logarithmic scale and on a decimal scale
Eb_N0_dB = 0 : 20;
Eb_N0 = 10.^(Eb_N0_dB/10);
%  Power of noise
N0 = Eb ./ Eb_N0;
% Define error numbers (used to calculate BER)
errors = zeros(1, length(Eb_N0_dB));  
% Define BER and theoretical BER
BER = zeros(1, length(Eb_N0_dB)); 
BER_Theoretical = zeros(1, length(Eb_N0_dB));

% Define generator matrix G and parity-check matrix H
G = [1, 1, 0, 1, 0, 0, 0;
        0, 1, 1, 0, 1, 0, 0;
        1, 1, 1, 0, 0, 1, 0;
        1, 0, 1, 0, 0, 0, 1];
    
H = [ 1, 0, 0, 1, 0, 1, 1;
         0, 1, 0, 1, 1, 1, 0;
         0, 0, 1, 0, 1, 1, 1];

%% Encoding and modulation
% The process of encoding is quite straightforward. We can directly
% multiply signal by the generator matrix as described in the assignment document (because Hamming codes are linear codes).
for i = 0 : data_length / k - 1
    data_encoded(1, i*7+1) = xor(data(1, i*4+1), xor(data(1, i*4+3), data(1, i*4+4)));    
    data_encoded(1, i*7+2) = xor(data(1, i*4+1), xor(data(1, i*4+2), data(1, i*4+3)));
    data_encoded(1, i*7+3) = xor(data(1, i*4+2), xor(data(1, i*4+3), data(1, i*4+4)));
    data_encoded(1, i*7+4) = data(1, i*4+1);
    data_encoded(1, i*7+5) = data(1, i*4+2);
    data_encoded(1, i*7+6) = data(1, i*4+3);
    data_encoded(1, i*7+7) = data(1, i*4+4);
end
% In order to do modulation with BPSK, we will play a tiny trick to let 0 be -1 and 1 be 1
data_mod = (data_encoded - 1/2) * 2;

%% Demodulation and decoding
for j = 1 : length(Eb_N0_dB)
%     Add interference of AWGN to the signal
    noise = sqrt(N0(j)/2) * randn(1, encoded_data_length);
    data_received = data_mod + noise;
%     Define demodulated signal and decoded signal
    data_demod = zeros(1, encoded_data_length);
    data_decoded = zeros(1, encoded_data_length);
    
%     Demodulate received signal: positive to 1, negative value to 0
    for d = 1 : encoded_data_length
        if (data_received(d) >= 0)
            data_demod(d) =1;
        else
            data_demod(d) = 0;
        end
    end
    
%   Syndrome
   s = zeros(1, 3);
    for t = 0 : data_length/k - 1
        s(1) = xor(data_demod(1,t*7+1),xor(data_demod(1,t*7+4),xor(data_demod(1,t*7+6),data_demod(1,t*7+7)))); 
        s(2) = xor(data_demod(1,t*7+2),xor(data_demod(1,t*7+4),xor(data_demod(1,t*7+5),data_demod(1,t*7+6))));
        s(3) = xor(data_demod(1,t*7+3),xor(data_demod(1,t*7+5),xor(data_demod(1,t*7+6),data_demod(1,t*7+7))));
        
        Correct signal with Decoding table for the (7, 4) Hamming code
        switch num2str(s)
            case num2str([0, 0, 1])
                data_demod(1, t*7+3) = ~data_demod(1, t*7+3);
            case num2str([0, 1, 0])
                data_demod(1, t*7+2) = ~data_demod(1, t*7+2);
            case num2str([0, 1, 1])
                data_demod(1, t*7+5) = ~data_demod(1, t*7+5);
            case num2str([1, 0, 0])
                data_demod(1, t*7+1) = ~data_demod(1, t*7+1);
            case num2str([1, 0, 1])
                data_demod(1, t*7+7) = ~data_demod(1, t*7+7);
            case num2str([1, 1, 0])
                data_demod(1, t*7+4) = ~data_demod(1, t*7+4);
            case num2str([1, 1, 1])
                data_demod(1, t*7+6) = ~data_demod(1, t*7+6);
        end
        
%         if (s == [0, 0, 1])
%             data_demod(1, t*7+3) = ~data_demod(1, t*7+3);
%         elseif (s == [0, 1, 0])
%                 data_demod(1, t*7+2) = ~data_demod(1, t*7+2);
%         elseif (s == [0, 1, 1])
%                 data_demod(1, t*7+5) = ~data_demod(1, t*7+5);
%         elseif (s == [1, 0, 0])
%                 data_demod(1, t*7+1) = ~data_demod(1, t*7+1);
%         elseif (s == [1, 0, 1])
%                 data_demod(1, t*7+7) = ~data_demod(1, t*7+7);
%         elseif (s == [1, 1, 0])
%                 data_demod(1, t*7+4) = ~data_demod(1, t*7+4);
%         elseif (s == [1, 1, 1])
%                 data_demod(1, t*7+6) = ~data_demod(1, t*7+6);
%         end
        
        
%         Decode the demodulated signal
        data_decoded(1,t*4+1) = data_demod(1,t*7+4);
        data_decoded(1,t*4+2) = data_demod(1,t*7+5);
        data_decoded(1,t*4+3) = data_demod(1,t*7+6);
        data_decoded(1,t*4+4) = data_demod(1,t*7+7);
    end
%     Count error number
    for b = 1 : data_length
        if(data_decoded(b) ~= data(b))
            errors(j) = errors(j) + 1;
        end
    end
%     Actural BER calculation
    BER(j) = errors(j) ./ data_length;
%     Theretical value
    BER_Theoretical(j) = qfunc(sqrt(2*Eb_N0(j)));
end

%% Plotting 
figure(1)
semilogy(Eb_N0_dB, BER, 'M-X', Eb_N0_dB, BER_Theoretical, 'B-O');                                       
xlabel('Eb/N0 (dB)');                           
ylabel('BER');
axis([0 20 10^-7 1])     
title('Bit Error Rate (BER) - Eb/N0');
legend('BER','Theoretial BER (Uncoded)');
grid on;  

%% NOTE
% The code is slow, but I have no time to profile the code. I guess the
% process of vector2str conversion takes too much time. We shall use
% if-else if-else to solve this. But in this code we can simply decrease
% the number of bit on input signal.











