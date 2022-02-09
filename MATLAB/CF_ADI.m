function [Z,k] = CF_ADI(A,C1,N,itmax,tol)
%   [Z,K] = CF_ADI(A,C1,N,ITMAX,TOL) risolve l'equazione di Lyapunov
%   A*X + X*A' + C1*C1' = 0 per mezzo dell'algoritmo Cholesky Factor ADI.
%   INPUT:
%       - A,C1: matrici dei coefficienti
%       - N: numero dei parametri di shift da adoperare
%       - ITMAX: numero massimo di iterazioni
%       - TOL: valore soglia per la condizione di arresto
%   OUTPUT:
%       - Z: matrice tale che X = Z*Z' risolve l'equazione
%       - K: numero di iterazioni svolte

[n,~] = size(A);
normA = norm(A,'fro');

%Iterazione 1
s = ADI_Suboptimal(A,N,2*N,2*N);
Z = sqrt(-2*real(s(1)))*((A+s(1)*eye(n))\C1);
U = Z;
normZ = norm(Z,'fro');

%Iterazioni successive alla prima
for k = 2:itmax
    m = mod(k,N);
    if (m == 0)
        m = N;
    end
    prec = m-1;
    if (prec == 0)
        prec = N;
    end
    U = sqrt(-2*real(s(m)))/sqrt(-2*real(s(prec)))*(U-(s(m)+conj(s(prec)))*((A+s(m)*eye(n))\U));
    Z = [Z,U];
    %Aggiornamento della norma della soluzione e valutazione della condizione di arresto
    normU = norm(U,'fro');
    normZ = sqrt(normZ^2 + normU^2);
    if (normU/(normZ*normA*2) < tol)
        return
    end
end
disp('Warning: raggiunto il numero massimo di iterazioni');