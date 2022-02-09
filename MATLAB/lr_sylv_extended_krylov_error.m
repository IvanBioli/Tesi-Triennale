function [X,k,error1,error2] = lr_sylv_extended_krylov_error(A,B,C1,C2,itmax,tol,flag)
%   [X,K,ERROR1,ERROR2] = LR_SYLV_EXTENDED_KRYLOV_ERROR(A,B,C1,C2,ITMAX,TOL,FLAG) risolve
%   l'equazione di Sylvester A*X + X*B = C1*C2' per mezzo dell'algoritmo 
%   proiettivo basato sugli Spazi di Krylov Estesi. Funzione finalizzata 
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
%         indicato nel paper di Heyouni

%Operazioni preliminari
[~,s] = size(C1);
[LA,UA] = lu(A);
[LB,UB] = lu(B');
if  (flag ~= 0)
    C = C1*C2';
else
    error1 = NaN;
end
H = zeros(2*s);
itmin = 20;
normA = norm(A,'fro');
normB = norm(B,'fro');

%Calcolo di V1 e W1
[V,~] = qr([A*C1,UA\(LA\C1)],0);
[W,~] = qr([B'*C2,UB\(LB\C2)],0);

for k = 1:itmax
    %Formulazione e risoluzione del problema di taglia ridotta
    Atilde = V'*A*V;
    Btilde = W'*B*W;
    C1tilde = V'*C1;
    C2tilde = W'*C2;
    Ctilde = C1tilde*C2tilde';
    Y = sylv_bartels_stewart(Atilde,Btilde,Ctilde);
    
    %Controllo della convergenza
    if  (flag ~= 0)
        X = V*Y*W';
        error1(k) = norm(A*X+X*B-C,'fro')/(norm(Y,'fro')*(normA+normB));
    end
    m = k-1;
    if k>1
        error2(m) = sqrt(norm(Atilde(2*m*s+1:2*m*s+s,2*(m-1)*s+1:2*m*s)*Y(2*(m-1)*s+1:2*m*s,:),'fro')^2+norm(Btilde(2*(m-1)*s+1:2*m*s,2*m*s+1:2*m*s+s)'*Y(:,2*(m-1)*s+1:2*m*s)','fro')^2)/(norm(Y,'fro')*(normA+normB));
        if (error2(k-1) < tol)% ||(m>1 && error2(m-1)<error2(m) && k>itmin))
            X = V*Y*W';
            return
        end
    end
    
    %Aggiornamento dello spazio di Krylov
    %Partizionamento di V e W per ottenere Vtilde e Wtilde
    Vplus = V(:,end-2*s+1:end-s);
    Vmin = V(:,end-s+1:end);
    Vtilde = [A*Vplus,UA\(LA\Vmin)];
    Wplus = W(:,end-2*s+1:end-s);
    Wmin = W(:,end-s+1:end);
    Wtilde = [B'*Wplus,UB\(LB\Wmin)];
    
    %Ortogonalizzazione di Vtilde risepetto a V e di Wtilde rispetto a W
    %Poichè le colonne di V e W sono ortonormali conviene calcolare come
    %sotto invece che con:
    %   Vtilde = Vtilde-V*V'*Vtilde
    %   Wtilde = Wtilde-W*W'*Wtilde
    for i = 1:k
        H = V(:,((i-1)*2*s+1):(i*2*s))'*Vtilde;
        Vtilde = Vtilde - V(:,((i-1)*2*s+1):(i*2*s))*H;
        H = W(:,((i-1)*2*s+1):(i*2*s))'*Wtilde;
        Wtilde = Wtilde - W(:,((i-1)*2*s+1):(i*2*s))*H;
    end
    
    %Aggiornamento di V e W
    [Vtilde,~]=qr(Vtilde,0);
    V = [V,Vtilde];
    [Wtilde,~]=qr(Wtilde,0);
    W = [W,Wtilde];
end
disp('Warning: raggiunto il numero massimo di iterazioni');
X = V(:,1:end-2*s)*Y*W(:,1:end-2*s)';
