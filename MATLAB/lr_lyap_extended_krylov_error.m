function [X,k,error1,error2] = lr_lyap_extended_krylov_error(A,C1,itmax,tol,flag)
%   [X,K,FLAG] = LR_LYAP_EXTENDED_KRYLOV(A,C1,ITMAX,TOL,FLAG) risolve
%   l'equazione di Lyapunov A*X + X*A' = C1*C1' per mezzo dell'algoritmo 
%   proiettivo basato sugli Spazi di Krylov Estesi.
%   INPUT:
%       - A,B,C1: matrici dei coefficienti
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
if  (flag ~= 0)
    C = C1*C1';
else
    error1 = NaN;
end
H = zeros(2*s);
itmin = 20;
normA = norm(A,'fro');

%Calcolo di V1 e W1
[V,~] = qr([A*C1,UA\(LA\C1)],0);

for k = 1:itmax
    %Formulazione e risoluzione del problema di taglia ridotta
    Atilde = V'*A*V;
    C1tilde = V'*C1;
    Ctilde = C1tilde*C1tilde';
    Y = lyap_bartels_stewart(Atilde,Ctilde);
    
    %Controllo della convergenza
    if  (flag ~= 0)
        X = V*Y*V';
        error1(k) = norm(A*X+X*A'-C,'fro')/(norm(Y,'fro')*2*normA);
    end
    m = k-1;
    if k>1
        error2(m) = sqrt(2)*norm(Atilde(2*m*s+1:2*m*s+s,2*(m-1)*s+1:2*m*s)*Y(2*(m-1)*s+1:2*m*s,:),'fro')/(norm(Y,'fro')*2*normA);
        if (error2(k-1) < tol ||(m>1 && error2(m-1)<error2(m) && k>itmin))
            X = V*Y*V';
            return
        end
    end
    
    %Aggiornamento dello spazio di Krylov
    %Partizionamento di V per ottenere Vtilde
    Vplus = V(:,end-2*s+1:end-s);
    Vmin = V(:,end-s+1:end);
    Vtilde = [A*Vplus,UA\(LA\Vmin)];
    
    %Ortogonalizzazione di Vtilde risepetto a V
    %Poichè le colonne di V sono ortonormali conviene calcolare come
    %sotto invece che con:
    %   Vtilde = Vtilde-V*V'*Vtilde
    for i = 1:k
        H = V(:,((i-1)*2*s+1):(i*2*s))'*Vtilde;
        Vtilde = Vtilde - V(:,((i-1)*2*s+1):(i*2*s))*H;
    end
    
    %Aggiornamento di V
    [Vtilde,~]=qr(Vtilde,0);
    V = [V,Vtilde];
end
disp('Warning: raggiunto il numero massimo di iterazioni');
X = V(:,1:end-2*s)*Y*V(:,1:end-2*s)';
