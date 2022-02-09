function [X,k,error1,error2] = lr_sylv_krylov_error(A,B,C1,C2,itmax,tol,flag)
%   [X,K,ERROR1,ERROR2] = LR_SYLV_KRYLOV_ERROR(A,B,C1,C2,ITMAX,TOL,FLAG) risolve 
%   l'equazione di Sylvester A*X + B*X = C1*C2' per mezzo dell'algoritmo 
%   proiettivo basato sugli Spazi di Krylov Standard. Funzione finalizzata 
%   all'analisi della convergenza in termini di numero di iterazioni svolte
%   e precisione raggiunta, NON attendibile dal punto di vista del tempo di
%   esecuzione.
%   INPUT:
%       - A,B,C1,C2: matrici dei coefficienti
%       - ITMAX: numero massimo di iterazioni
%       - TOL: valore soglia per la condizione di arresto
%       - FLAG: se TRUE viene effettivamente calcolato ERROR1, altrimenti
%         ERROR1 viene impostato di default a NaN
%   OUTPUT:
%       - X: soluzione dell'equazione
%       - K: numero di iterazioni svolte
%       - ERROR1: vettore contenete l'errore relativo commesso sulla
%         approssimazione della soluzione ad ogni passo, cioè per ogni
%         passo norm(A*X+X*B-C,'fro')/norm(X,'fro'), dove X = V*Y*W'
%       - ERROR2: vettore contenete l'errore "calcolato velocemente" come
%         indicato nel paper di Simoncini

%Operazioni preliminari
[~,s] = size(C1);
[V,~] = qr(A*C1,0);
[W,~] = qr(B'*C2,0);
if  (flag ~= 0)
    C = C1*C2';
else
    error1 = NaN;
end
normA = norm(A,'fro');
normB = norm(B,'fro');
    
for k = 1:itmax
    %Formulazione e risoluzione del problema di taglia ridotta
    C1tilde = V'*C1;
    C2tilde = W'*C2;
    Ctilde = C1tilde*C2tilde';
    
    %Controllo della convergenza
    if (k > 1)
        error2(k-1) = norm(HV*Y*[eye(k-1),zeros(k-1,s)]+[eye(k-1);zeros(s,k-1)]*Y*HW'-Ctilde,'fro')/(norm(Y,'fro')*(normA+normB));
        if (error2(k-1) < tol)
            X = V(:,1:end-s)*Y*W(:,1:end-s)';
            return
        end
    end
    
    Atilde = V'*A*V;
    Btilde = W'*B*W;
    Y = sylv_bartels_stewart(Atilde,Btilde,Ctilde);
    if  (flag ~= 0)
        X = V*Y*W';
        error1(k) = norm(A*X+X*B-C,'fro')/(norm(X,'fro')*(normA+normB));
    end
    
    %Aggiornamento dello spazio di Krylov
    %Partizionamento di V e W per ottenere Vtilde e Wtilde
    Vtilde = A*V(:,end-s+1:end);
    Wtilde = B'*W(:,end-s+1:end);
    
    %Ortogonalizzazione di Vtilde risepetto a V e di Wtilde rispetto a W
    %Poichè le colonne di V e W sono ortonormali conviene calcolare come
    %sotto invece che con:
    %   Vtilde = Vtilde-V*V'*Vtilde
    %   Wtilde = Wtilde-W*W'*Wtilde
    for i = 1:k
        HV(((i-1)*s+1):(i*s),((k-1)*s+1):(k*s)) = V(:,((i-1)*s+1):(i*s))'*Vtilde;
        Vtilde = Vtilde - V(:,((i-1)*s+1):(i*s))*HV(((i-1)*s+1):(i*s),((k-1)*s+1):(k*s));
        HW(((i-1)*s+1):(i*s),((k-1)*s+1):(k*s)) = W(:,((i-1)*s+1):(i*s))'*Wtilde;
        Wtilde = Wtilde - W(:,((i-1)*s+1):(i*s))*HW(((i-1)*s+1):(i*s),((k-1)*s+1):(k*s));
    end
    
    %Aggiornamento di V e W
    [Vtilde,HV((k*s+1):((k+1)*s),((k-1)*s+1):(k*s))]=qr(Vtilde,0);
    V = [V,Vtilde];
    [Wtilde,HW((k*s+1):((k+1)*s),((k-1)*s+1):(k*s))]=qr(Wtilde,0);
    W = [W,Wtilde];
end
disp('Warning: raggiunto il numero massimo di iterazioni');
X = V(:,1:end-s)*Y*W(:,1:end-s)';
