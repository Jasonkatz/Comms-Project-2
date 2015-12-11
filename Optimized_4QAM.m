% A skeleton BER script for a wireless link simulation
% Illustrates how to simulate modulation and compare it to a theoretical
% BER
clear all;close all;clc
% For the final version of this project, you must use these 3
% parameter. You will likely want to set numIter to 1 while you debug your
% link, and then increase it to get an average BER.
totpak = 50     %totl number of packets sent
nSym = 1000;    % The number of symbols per packet
M = 4;        % The M-ary number, 2 corresponds to binary modulation
k=log2(M); %M=2^k
T = 4; %samples per symbol
train = 100; %training bits
SNR_Vec = 1:12;
EbNo = SNR_Vec - 10*log10(k/T); %calculate EbNo for 
%Es/N0 (dB)=10log10(0.5Tsym/Tsamp)+SNR? (dB) for real input signals
lenSNR = length(SNR_Vec);
h = waitbar(0,'Initializing waitbar...');

%chan = 1;          % No channel
chan = [1 .2 .4]; % Somewhat invertible channel impulse response, Moderate ISI
%chan = [0.227 0.460 0.688 0.460 0.227]';   % Not so invertible, severe ISI

% Create a vector to store the BER computed during each iteration
berVec = zeros(1, lenSNR);
BERVEC = zeros(1, lenSNR);

for packetnum = 1:totpak;
    
bits = randint(1, nSym*k, [0 1]);% Generate random bits 
% xx changed M to k
% New bits must be generated at every
% iteration

% If you increase the M-ary number, as you most likely will, you'll need to
% convert the bits to integers. See the BIN2DE function
% For binary, our MSG signal is simply the bits


msg = bi2de(reshape(bits,k,nSym).').'; %convert to base M ints
%msg = bits;


tx = rectpulse(qammod(msg,M,0,'gray'),T);% BPSK modulate the signal

    if isequal(chan,1)
        txChan = tx;
    else
        txChan = filter(upsample(chan,T),1,tx); % Apply the channel
        txChan(1) = txChan(1) +i*10^-7 ;
    end

for j = 1:lenSNR % one iteration of the simulation at each SNR Value      
    txNoisy = awgn(txChan,SNR_Vec(j),'measured'); % Add AWGN
    forgetfactor = .998;
    alg = rls(forgetfactor,.06);
    eqobj = dfe(17, 9, alg, qammod(0:M-1,M,0,'gray'));
    txeq = equalize(eqobj, txNoisy,tx(1:train)); % equalize using dfe
    %txeq = txNoisy;
    txeqid = myintdump(txeq,T,2,T);
    rx = qamdemod(txeqid,M,0,'gray'); % Demodulate
    %eyediagram(rx,T);
   
    % Again, if M was a larger number, I'd need to convert my symbols
    % back to bits here.
    rxMSG = reshape(de2bi(rx(train+1:end),k).',1,(nSym-train)*k);
    
    % Compute and store the BER for this iteration
    
    [biterrors berVec(j)] = biterr(bits(train*k+1:end), rxMSG);  % We're interested in the BER, which is the 2nd output of BITERR
   %waitbar(j/(lenSNR),h,sprintf('%d%% along...',100*j/lenSNR));
    waitbar(packetnum/totpak,h,sprintf('%d%% along...',100*packetnum/totpak));
BERVEC(j) = berVec(j) + BERVEC(j);
end
end
berVec = BERVEC/packetnum;
sPlotFig = scatterplot(txeqid(train+1:end),1,0,'g.');
hold on;
scatterplot(qammod(0:M-1,M,0,'gray'),1,0,'k*',sPlotFig)
nSym
packetnum
bitrate = ((nSym-train)*k-biterrors)/nSym
berVec
if lenSNR ~= 1
    figure;
    semilogy(SNR_Vec, berVec,'*')

    % Compute the theoretical BER for this scenario
    berTheory = berawgn(EbNo,'qam',M,'nondiff');
    hold on
    semilogy(SNR_Vec,berTheory,'r')
    xlabel('SNR (dB)')
    ylabel('BER (errbits/informationbits)')
    title('4-QAM SNR to BER')
    legend('BER', 'Theoretical BER','Location','best')
end

close(h);