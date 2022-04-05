clc
clear 
close all

% Christian Kamps - 101004691

% INCLUDED FILES:
%   - root_raised_cosine.m : A given RRC signal function
%   - channel.p : The given channel

% DO NOT RUN EVERYTHING AT THE SAME TIME:
%   - Certain components are disabled for simulation speed or have conditional requirements.
%   - In general, reduce the number of frames (Nf) to increase simulation speed   
%   - The current simulation setup is for various BER values for a RRC pulse

% If you are looking at PSD, disable BER calculations:
%   - Disable SNR loop
%   - Disable Channel
%   - Disable Reciever 
%   - Choose whether you are looking at multiple PSD or singular (NOT BOTH)
%   - Disable BER plots

% If you are looking at BER, disable PSD calculations:
%   - Disable all PSD calculations
%   - Disable PSD plot

% Loop for variable parameters:
%   - This is only for the figures generated for the report.
%   - Several plots at once for varying parameter values.
%   - These simulations take a while due to repetitive simulations loops
%   - Disable if not looking at varying paramter values.

%%%%%%%%%%%%%%%%%%%%%%
% System Parameters
%%%%%%%%%%%%%%%%%%%%%%
kmax = 500;                         % Maximum Delay
Nf = 1000;                          % Number of frames
Na = 1200;                          % Message length (bits)
T = 0.0035;                         % Symbol duration (seconds/symbol)
eta = 64;                           % Oversampling ratio (samples/symbol)
fc = 850;                           % Carrier frequency (Hz)
Ts = T / eta;                       % Sampling period (seconds/sample)
SNRdB = 0:10;	                    % SNR (Eb/No) range of interest (in dB)
totErrs = zeros(1, length(SNRdB));  % Total errors vector allocation

%%%%%%%%%%%%%%%%%%%%%%
% Symbol Mapping
%%%%%%%%%%%%%%%%%%%%%%
SM = [ +1+j +1-j -1+j -1-j];    % 4-QAM symbol map
M = length(SM);                 % Constellation size
bitsPerSymbol = log2(M);        % Number of message bits per symbol
Es = sum(abs(SM).^2) / M;       % Energy per symbol
Eb = Es / bitsPerSymbol;        % Energy per bit

%%%%%%%%%%%%%%%%%%%%%%
% Variable Parameters
%%%%%%%%%%%%%%%%%%%%%%
% L_v = [1 4 8 12 16];                % Pulse length values
% B_v = [0 0.2 0.4 0.6 0.8 1];        % Roll-off factor values
% T_v = [0.01 0.005 0.001 0.0005];    % Symbol period values
% Fc_v = [650 750 850 950 1050];      % Frequency Value
% Ideal_T_v = 0.0030:0.0001:0.0040;   % Ideal range of Symbol period values
% param = Symbol_period_vector;       % Select vector for loop

%%%%%%%%%%%%%%%%%%%%%%
% Loop for Variable Paramters
%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:length(param)
%     PSD_vector = [];
%     fc = param(i);
%     L = param(i)    
%     beta = param(i)
%     T = param(i); % If you enable variable T values you must enable Ts re-calculation
%     Ts = T / eta;  

%%%%%%%%%%%%%%%%%%%%%%
% Rectangular Pulse Shape
%%%%%%%%%%%%%%%%%%%%%%
% L = 1;                            % Must set L=1
% hT = sqrt(1/T) * ones(1, eta);    % Rectangular Pulse Shape 

%%%%%%%%%%%%%%%%%%%%%%
% RRC Pulse Shape
%%%%%%%%%%%%%%%%%%%%%%
L = 16;                                     % Truncated Pulse Duration (symbol periods)
beta = 1;                                   % Roll-off factor
hT = root_raised_cosine(beta, L, T, eta);   % Root-Raised Cosine Pulse Shape

%%%%%%%%%%%%%%%%%%%%%%
% Matched Filter
%%%%%%%%%%%%%%%%%%%%%%
hR = fliplr(hT);    % Receiver filter impulse response

%%%%%%%%%%%%%%%%%%%%%%
% Training Sequence
%%%%%%%%%%%%%%%%%%%%%%
Nt = 100;                               % Training Sequence Length
v_train = 1 - 2*randi([0, 1], 1, Nt);   % Training Symbols vector

%%%%%%%%%%%%%%%%%%%%%%
% Loop for SNR Values
%%%%%%%%%%%%%%%%%%%%%%
for nn = 1:length(SNRdB)
     nErrs = 0;                         % Reset total errors of a single message
     No = Eb * 10.^(-SNRdB(nn)/10);     % Single-sided noise PSD
   
    %%%%%%%%%%%%%%%%%%%%%%
    % Loop for Messages
    %%%%%%%%%%%%%%%%%%%%%%
    for nf = 1:Nf

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                     %
        % Transmitter                                                         %
        %                                                                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%
            % Message source 
            %%%%%%%%%%%%%%%%%%%%%%   
            a = randi([0 1], 1, Na);  % Generate message of bits of Na length
    
            %%%%%%%%%%%%%%%%%%%%%%
            % Symbol Mapper 
            %%%%%%%%%%%%%%%%%%%%%%            
            c = reshape(a, bitsPerSymbol, []);      % Group message bits into groups of 'bitsPerSymbol' bits
            ci = 2.^(bitsPerSymbol-1:-1:0) * c;     % Convert each group of bits into an int
            v = SM(ci+1);                           % Symbol map the groups of bits
            v_p = [v_train v];                      % Add training symbols to transmitted symbols
            
            %%%%%%%%%%%%%%%%%%%%%%
            % Pulse Shaping Filter 
            %%%%%%%%%%%%%%%%%%%%%%  
            vt = conv(upsample(v_p, eta), hT);  % Transmit Filter
            vt = vt(1:end-eta+1);               % Truncate extra samples
    
            %%%%%%%%%%%%%%%%%%%%%%
            % Modulator 
            %%%%%%%%%%%%%%%%%%%%%% 
            t = (0:length(vt)-1)*Ts;                        % Sample times
            vct = real(vt .* sqrt(2) .* exp(j*2*pi*fc*t));  % Modulated Signal       

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                     %
        % Channel                                                             %
        %                                                                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%
            % Bandpass AWGN Channel
            %%%%%%%%%%%%%%%%%%%%%%  
            wct = sqrt(1/Ts*No/2)*randn(1, length(vct)+kmax);   % Noise Signal
            rct = channel(vct, Ts) + wct;                       % Received Signal

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                     %
        % Receiver                                                            %
        %                                                                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%
            % Demodulator
            %%%%%%%%%%%%%%%%%%%%%%  
            t = (0:length(rct)-1)*Ts;                   % Sample times
            rot = rct .* sqrt(2) .* exp(-j*2*pi*fc*t);  % Demodulated Signal
    
            %%%%%%%%%%%%%%%%%%%%%%
            % Matched Filter
            %%%%%%%%%%%%%%%%%%%%%%  
            rt = conv(rot, hR) * Ts;    % Matched filter signal
    
            %%%%%%%%%%%%%%%%%%%%%%
            % Clock Recovery
            %%%%%%%%%%%%%%%%%%%%%%          
            k_vector = 0:500 ;                                                      % k-vector representing possible delays           
            trainining_symbol_index = (L*eta:eta:(L*eta+eta*(Nt-1))) + k_vector';   % All possible training symbol indexs for any delay (from 0 to 500): L𝜂+n𝜂+k      
            r_oversample = rt(trainining_symbol_index);                             % Oversample of the matched filter output r(t): 𝑟̃_𝑛
            mu_k = (1/Nt).*sum( r_oversample.*v_train , 2);                         % Correlation coefficient for the delay
            [max_value,index_value] = max(mu_k);                                    % Max correlation coefficient value and its k value
            k = index_value-1;                                                      % The k value (propagation delay) for the max correlation coefficient
            
            %%%%%%%%%%%%%%%%%%%%%%
            % Carrier Recovery
            %%%%%%%%%%%%%%%%%%%%%%                  
            r = rt(k+L*eta:eta:end);                % Sampler
            r_train = r(1:Nt);                      % Training sequence samples
            r_cut = r(Nt+1:Nt+Na/bitsPerSymbol);    % Discard training sequence samples          
            q = (1/Nt) * sum(r_train./v_train);     % Calculate q for phase error
            z = r_cut / q;                          % Compensate for phase error 
    
            %%%%%%%%%%%%%%%%%%%%%%
            % 4-QAM Decision Device
            %%%%%%%%%%%%%%%%%%%%%%  
            ah = zeros(1, Na);              % Pre-allocate received message
            ah(1:2:Na) = (real(z) < 0);     % If the symbol real component is less than 0 the second bit is 1
            ah(2:2:Na) = (imag(z) < 0);     % If the symbol imag component is less than 0 the first bit is 1
    
            %%%%%%%%%%%%%%%%%%%%%%
            % Data Sink
            %%%%%%%%%%%%%%%%%%%%%%  
            nErrs = nErrs + sum(xor(a, ah));    % Count errors  
            
        %%%%%%%%%%%%%%%%%%%%%%
        % PSD
        %%%%%%%%%%%%%%%%%%%%%%    
%         Ns = length(vct);               % Length of bandpass signal
%         Vcf = fftshift(fft(vct));       % Calculate FFT
%         PSD = Vcf.*conj(Vcf) * Ts / Ns; % Calculate PSD
%         f = (-Ns/2:Ns/2-1) / (Ns*Ts);   % Frequencies        
%         PSD_vector(nf,:) = PSD;         % Store every PSD value in a matrix
    end
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Bit Error Rate (BER) for multiple SNR values
    %%%%%%%%%%%%%%%%%%%%%%   
    totErrs(nn) = nErrs;            % Total errors vector
    BER_sim = totErrs / (Na*Nf)     % Calculate the BER - Keep this as a print to view simulation progress
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Bit Error Rate (BER) for one SNR value (9)
    %%%%%%%%%%%%%%%%%%%%%%       
%     totErrs = nErrs;
%     BER_sim(:,i) = totErrs / (Na*Nf)

    %%%%%%%%%%%%%%%%%%%%%%
    % PSD with variable values
    %%%%%%%%%%%%%%%%%%%%%%   
    % Only enable this if the variable paramters for loop is enabled.
%     f_v{i,:} = f;                                                       % Frequency matrix
%     PSD_db{i,:} = 10*log10(sum(PSD_vector(:,(1:length(PSD))) )./Nf);    % Sum all PSD's and divide by the message amount to get an average
   
end

%%%%%%%%%%%%%%%%%%%%%%
% PSD without variable values
%%%%%%%%%%%%%%%%%%%%%%   
% PSD_db = 10*log10(sum(PSD_vector(:,(1:length(PSD))) )./Nf);                  % Sum all PSD's and divide by the message amount to get an average

%%%%%%%%%%%%%%%%%%%%%%
% Theoretical Bit Error Rate (BER) 
%%%%%%%%%%%%%%%%%%%%%%   
BER_theory = qfunc(sqrt(2*10.^(SNRdB/10))); % Theoretical BER for 4-QAM
% SNR = 10.^(SNRdB./10);                    % Power-ratio of SNR
% BER_theory = (3/8)*erfc(sqrt((2/5)*SNR))+(2/8)*erfc(3*sqrt((2/5)*SNR))-(1/8)*erfc(5*sqrt((2/5)*SNR)); % BER of 16-QAM

%%%%%%%%%%%%%%%%%%%%%%
% BER Plot
%%%%%%%%%%%%%%%%%%%%%%
figure(1)
semilogy(SNRdB, BER_theory, SNRdB, BER_sim,'LineWidth',2)
xlabel('E_b / N_0 (dB)');
ylabel('BER');
title('Bit Error Rate vs. SNR')
grid on
legend('Theoretical','Simulated');

% %%%%%%%%%%%%%%%%%%%%%%
% % PSD Plot
% %%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% hold on
% xlabel('Frequency (Hz)');
% ylabel('PSD(dB)');
% title('Average PSD of Designed Bandpass Signal v_c(t)')
% %axis([450 1250 -80 10]); % Main frequency component
% ylim([-80 10])
% xlim([-3000 3000])
% %legend('Rectangular Pulse', 'RRC Pulse');
% %legend('f_c=650', 'f_c=750','f_c=850','f_c=950','f_c=1050');
% %legend('L=1', 'L=4','L=8','L=12','L=16');
% %legend('Beta=0', 'Beta=0.2','Beta=0.4','Beta=0.6','Beta=0.8','Beta=1');
% %legend('T=0.01', 'T=0.005','T=0.001','T=0.0005');
% 
% %%%%%%%%%%%%%%%%%%%%%%
% % PSD Plot for no variable values
% %%%%%%%%%%%%%%%%%%%%%%    
% plot(f, PSD_db);
% 
% %%%%%%%%%%%%%%%%%%%%%%
% % PSD Plot for variable values
% %%%%%%%%%%%%%%%%%%%%%%    
% for i = 1:length(param)
%     plot(f_v{i,:}, PSD_db{i,:},'linewidth',2);
% end
% hold off