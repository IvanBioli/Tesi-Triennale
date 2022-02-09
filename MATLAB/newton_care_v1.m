function [X,k] = newton_care_v1(A,B,C,X0,tol,itmax)
%   [X,K] = NEWTON_CARE_V1(A,B,C,X0,TOL,ITMAX) risolve la CARE 
%   XA + A'X - XBX = C per mezzo del metodo di Newton basato sul calcolo di
%   H_K = X_(K+1) - X_K
%   INPUT
%       - A, B, C: matrici dei coefficienti
%       - X0: approssimazione iniziale
%       - TOL: tolleranza per la condizione di arresto
%       - ITMAX: numero massimo di iterazioni
%   OUTPUT
%       - X : soluzione della CARE
%       - K: numero di iterazioni del metodo di Newton effettuate

X = X0;
err = 1;
k = 0;
normA = norm(A,'fro');
normB = norm(B,'fro');

while (err > tol && k < itmax)
    RX = X*A + A'*X - X*B*X - C;
    H = lyap_bartels_stewart(A'-X*B,-RX);
    X = X + H;
    err = norm(H,'fro')/(norm(X,'fro')*(normA+normB));
    k = k + 1;
end
if k == itmax
    disp('Warning: raggiunto il numero massimo di iterazioni')
end

