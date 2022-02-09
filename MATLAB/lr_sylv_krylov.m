function [X,k] = lr_sylv_krylov(A,B,C1,C2,itmax,tol)
%   [X,K] = LR_SYLV_KRYLOV(A,B,C1,C2,ITMAX,TOL) risolve 
%   l'equazione di Sylvester A*X + B*X = C1*C2' per mezzo dell'algoritmo 
%   proiettivo basato sugli Spazi di Krylov Standard.
%   INPUT:
%       - A,B,C1,C2: matrici dei coefficienti
%       - ITMAX: numero massimo di iterazioni
%       - TOL: valore soglia per la condizione di arresto
%   OUTPUT:
%       - X: soluzione dell'equazione
%       - K: numero di iterazioni svolte

%Operazioni preliminari
[~,s] = size(C1);
[V,~] = qr(A*C1,0);
[W,~] = qr(B'*C2,0);
normA = norm(A,'fro');
normB = norm(B,'fro');

for k = 1:itmax
    %Formulazione e risoluzione del problema di taglia ridotta
    C1tilde = V'*C1;
    C2tilde = W'*C2;
    Ctilde = C1tilde*C2tilde';
    %Controllo della convergenza
    if (k > 1)
        err = norm(HV*Y*[eye((k-1)*s),zeros((k-1)*s,s)]+[eye((k-1)*s);zeros(s,(k-1)*s)]*Y*HW'-Ctilde,'fro')/(norm(Y,'fro')*(normA+normB));
        if (err < tol)
            X = V(:,1:end-s)*Y*W(:,1:end-s)';
            return
        end
    end
    
    Atilde = V'*A*V;
    Btilde = W'*B*W;
    Y = sylv_bartels_stewart(Atilde,Btilde,Ctilde);
    
    %Aggiornamento dello spazio di Krylov
    %Partizionamento di V e W per ottenere Vtilde e Wtilde
    Vtilde = A*V(:,end-s+1:end);
    Wtilde = B'*W(:,end-s+1:end);
    
    %Ortogonalizzazione di Vtilde risepetto a V e di Wtilde rispetto a W
    %Poiche' le colonne di V e W sono ortonormali conviene calcolare come
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
