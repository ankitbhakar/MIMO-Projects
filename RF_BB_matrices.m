function[RFmat, BBmat] = RF_BB_matrices(numAnt,numRF,N_Beam)
NBlk = numAnt/numRF;
RFmat = 1/sqrt(numAnt)*dftmtx(numAnt);
U = random_unitary(numRF);
V = random_unitary(N_Beam/NBlk);
BB_diag = U*([eye(N_Beam/NBlk) zeros(N_Beam/NBlk,numRF-N_Beam/NBlk)]')*V';
BBmat = kron(eye(NBlk),BB_diag); 
 