function [Z,k,error1,error2] = CF_ADI_error(A,C1,N,itmax,tol,flag)
%   [Z,K,ERROR1,ERROR2] = CF_ADI(A,C1,N,ITMAX,TOL,FLAG) risolve 
%   l'equazione di Lyapunov A*X + X*A' + C1*C1' = 0 per mezzo 
%   dell'algoritmo Cholesky Factor ADI. Funzione finalizzata all'analisi 
%   della convergenza in termini di numero di iterazioni svolte e 
%   precisione raggiunta, NON attendibile dal punto di vista del tempo di 
%   esecuzione.
%   INPUT:
%       - A,C1: matrici dei coefficienti
%       - N: numero dei parametri di shift da adoperare
%       - ITMAX: numero massimo di iterazioni
%       - TOL: valore soglia per la condizione di arresto
%       - FLAG: se TRUE viene effettivamente calcolato ERROR1, altrimenti
%         ERROR1 viene impostato di default a NaN
%   OUTPUT:
%       - Z: matrice tale che X = Z*Z' risolve l'equazione
%       - K: numero di iterazioni svolte
%       - ERROR1: vettore contenete l'errore relativo commesso sulla
%         approssimazione della soluzione ad ogni passo, cioè per ogni
%         passo norm(A*X+X*A'+C1*C1','fro')/norm(Z,'fro') dove X = Z*Z'
%       - ERROR2: vettore contenete l'errore "calcolato velocemente", cioè
%         per ogni passo norm(U,'fro')/norm(Z,'fro')

[n,~] = size(A);
if (flag == 0)
    error1 = NaN;
else
    C = C1*C1';
end
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
    X = Z*Z';
    
    %Calcolo di error1(k) e error2(k) e valutazione della condizione di arresto
    if (flag ~= 0)
        X = Z*Z';
        error1(k-1) = norm(A*X+X*A'+C,'fro')/(norm(X,'fro')*2*normA);
    end
    normU = norm(U,'fro');
    normZ = sqrt(normZ^2 + normU^2);
    error2(k-1) = normU/(normZ*2*normA);
    if (normU/(normZ*2*normA) < tol)
        return
    end
end
disp('Warning: raggiunto il numero massimo di iterazioni');