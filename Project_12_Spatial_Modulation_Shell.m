close all;
clear all;

SNRdB = 0:2:10;
ITER = 100000;
Nt = 4;
Nr = 4;
M = 2;
bpcu = log2(Nt*M);
BERopt = zeros(1, length(SNRdB));

for ite = 1:ITER
  TxBits = randi([0,1], 1, bpcu)
  H = (1/sqrt(2))*(randn(Nr,Nt) + 1j*randn(Nr,Nt))
  RxNoise = (1/sqrt(2))*(randn(Nr,1) + 1j*randn(Nr,1))
  iBit = (2*TxBits(1,bpcu)-1) % getting the last bit which represents BPSK
  antIndex = 1+bin2dec(num2str(TxBits(1,1:log2(Nt))))  % It picks out the first 2 bits to pick out the antenna index
  for K = 1:length(SNRdB)
    rho = 10^(SNRdB(K)/10)
    RxVec = sqrt(rho)*H(:,antIndex)*iBit + RxNoise

    % Optimal ML Decoder
    MLobj = sum(abs(sqrt(rho)*[-H,H] - repmat(RxVec,1,2*Nt)).^2)   % We are calc the distance between receive vector and all possible symbols
    [minVal,minIdx] = min(MLobj)   % Getting the index of the min distance
    BitDec = (minIdx > M*Nt/2)    # decoding the info bit (MNt/2 = 4)
    antIndDec = dec2bin(minIdx-1-BitDec*Nt,log2(Nt))~='0'
    DecBits = [antIndDec, BitDec]
    BERopt(K) = BERopt(K) + sum(DecBits ~= TxBits)
  endfor
endfor

BERopt = BERopt/(bpcu*ITER);

semilogy(SNRdB,BERopt,'b- s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',9.0)
hold on;
axis tight;
xlabel('SNR (dB)')
ylabel('P_e')
legend('BER','Location','SouthWest')
title('BER Performance of Spatial Modulation vs SNR')
grid on;
axis tight;
