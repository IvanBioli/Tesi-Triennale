function X = sylv_vectorial(A,B,C)
% X=SYLV_VECTORIAL(A,B,C) risolve l'equazione di Sylvester AX + XB = C sfruttando la formulazione come sistema lineare tramite prodotto di Kroenecker
% A, B, C: matrici dei coefficienti (di taglie compatibili)
% X soluzione
[m,n] = size(C);
H = kron(eye(n), A) + kron(B.', eye(m));
Cvec = reshape(C,m*n,1);
Xvec = H\Cvec;
X = reshape(Xvec,m,n);