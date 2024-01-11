function U = random_unitary(m)
H = (randn(m,m)+1j*randn(m,m))/sqrt(2);
[U,S,V] = svd(H);