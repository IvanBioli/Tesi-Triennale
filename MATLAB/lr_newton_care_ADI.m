function [X,k] = lr_newton_care_ADI(A,BU,C,X0,tol,itmax)
%   [X,K] = LR_NEWTON_CARE_ADI(A,BU,C,X0,TOL,ITMAX) risolve la CARE
%   X*A + A'*X - X*B*X = C per mezzo del metodo di Newton sfruttando la
%   fattorizzazione di rango basso di B = BU*BU' e il solutore di rango
%   basso basato sul metodo ADI.
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
normA = norm(A,'fro');
normB = norm(BU,'fro')^2;
Ak = A'-X0*BU*(BU');
Ck = C-X0*BU*(BU')*X0;
X = lyap_bartels_stewart(Ak,Ck);
H = X-X0;
N = 20;                        
itmaxADI = 250;
while (err > tol && k < itmax)
    Ak = A'-X*BU*(BU');
    Z = CF_ADI(Ak,H*BU,N,itmaxADI,tol);                           %Si usa il metodo ADI per l'update
    H = -Z*Z';
    X = X + H;
    err = norm(H,'fro')/(norm(X,'fro')*(normA+normB));
    k = k + 1;
end
if k == itmax
    disp('Warning: raggiunto il numero massimo di iterazioni')
end
