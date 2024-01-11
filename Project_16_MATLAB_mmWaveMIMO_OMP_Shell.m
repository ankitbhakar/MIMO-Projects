close all; 
clear all; 
rng('shuffle');

t = 32; r = 32; 
numRF = 8;  
N_Beam = 24; 
G = 32; 
ITER = 10; 
L = 5; 

omp_thrld = 10;
SNRdB = 10:10:50;
SNR = 10.^(SNRdB/10);
mseOMP = zeros(length(SNRdB),1);
mseGenie = zeros(length(SNRdB),1);

A_T = 1/sqrt(t)*exp(-j*pi*[0:t-1]'*[2/G*[0:G-1] - 1]); 
A_R = A_T;

[FRF, FBB] = RF_BB_matrices(t,numRF,N_Beam);
[WRF, WBB] = RF_BB_matrices(r,numRF,N_Beam);
Qtil = kron((FBB.')*(FRF.'),(WBB')*(WRF'));


for ix = 1:ITER
% Write the lines here

    for ux = 1:length(SNRdB)
        % Write the lines here           
    end
end
mseOMP = mseOMP/ITER; 
mseGenie = mseGenie/ITER;

semilogy(SNRdB,mseOMP,'b s-','linewidth',3.0,'MarkerFaceColor','b','MarkerSize',9.0);
hold on;
semilogy(SNRdB,mseGenie,'m o-.','linewidth',3.0,'MarkerFaceColor','m','MarkerSize',9.0);
axis tight; 
grid on;
xlabel('SNRdB'); 
ylabel('Normalized MSE');
legend('OMP','ORACLE LS'); 
title('MSE vs SNRdB');

