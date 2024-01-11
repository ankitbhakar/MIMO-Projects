close all;
clear all;

SNRdB = 0:2:10;
ITER = 100000;
Nt = 4;
Nr = 4;
bpcu = log2(Nt);
BERopt = zeros(1, length(SNRdB));

for ite = 1:ITER
  TxBits = randi([0,1], 1, bpcu)
  H = (1/sqrt(2))*(randn(Nr,Nt) + 1j*randn(Nr,Nt))
  RxNoise = (1/sqrt(2))*(randn(Nr,1) + 1j*randn(Nr,1))
  antIndex = 1+bin2dec(num2str(TxBits))  % It picks out the first 2 bits to pick out the antenna index
  for K = 1:length(SNRdB)
    rho = 10^(SNRdB(K)/10)
    RxVec = sqrt(rho)*H(:,antIndex) + RxNoise

    % Optimal ML Decoder
    MLobj = sum(abs(sqrt(rho)*H - repmat(RxVec,1,Nt)).^2)   % We are calc the distance between receive vector and all possible symbols
    [minVal,minIdx] = min(MLobj)   % Getting the index of the min distance
    DecBits = dec2bin(minIdx-1,log2(Nt))~='0'
    BERopt(K) = BERopt(K) + sum(DecBits ~= TxBits)
  endfor
endfor

BERopt = BERopt/(bpcu*ITER);

semilogy(SNRdB,BERopt, 'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',9.0)
grid;
axis tight;
xlabel('SNR (dB)')
ylabel('P_e')
legend('BER Optimal','Location','SouthWest')
title('BER Performance of SSK vs SNR')

