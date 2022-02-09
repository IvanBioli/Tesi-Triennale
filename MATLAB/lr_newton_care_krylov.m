function [X,k] = lr_newton_care_krylov(A,BU,C,X0,tol,itmax)
%   [X,K] = LR_NEWTON_CARE_KRYLOV(A,BU,C,X0,TOL,ITMAX) risolve la CARE
%   X*A + A'*X - X*B*X = C per mezzo del metodo di Newton sfruttando la
%   fattorizzazione di rango basso di B = BU*BU' e il solutore di rango
%   basso basato sugli Spazi di Krylov Estesi.
%   INPUT
%       - A, BU, C: matrici dei coefficienti
%       - X0: approssimazione iniziale
%       - tol: tolleranza per la condizione di arresto
%       - itmax: numero massimo di iterazioni
%   OUTPUT
%       - X : soluzione della CARE
%       - K: numero di iterazioni del metodo di Newton effettuate

err = 1;
k = 0;
Ak = A'-X0*BU*(BU');
Ck = C-X0*BU*(BU')*X0;
X = lyap_bartels_stewart(Ak,Ck);
H = X-X0;
itmaxkryl = 100;
while (err > tol && k < itmax)
    Ak = A'-X*BU*(BU');
    H = lr_lyap_extended_krylov(Ak,H*BU,itmaxkryl,tol);           %Si usa il metodo di Krylov per l'update
    % H = lyap_bartels_stewart(Ak,H*BU*(BU')*(H'));
    X = X + H;
    err = norm(H,1)/norm(X,1);
    k = k + 1;
end
if k == itmax
    disp('Warning: raggiunto il numero massimo di iterazioni')
end