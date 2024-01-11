function [alpha, Psi, A_R_genie, A_T_genie] = mmWaveMIMOChannelGenerator(A_R,A_T,G,L)

AoDlist = randperm(G,L); AoAlist = randperm(G,L);
A_T_genie = A_T(:,AoDlist); A_R_genie = A_R(:,AoAlist);
alpha = 1/sqrt(2)*(randn(L,1) + 1j*randn(L,1));
for I = 1:L
    Psi(:,I) = kron(conj(A_T(:,AoDlist(I))),A_R(:,AoAlist(I)));
end
