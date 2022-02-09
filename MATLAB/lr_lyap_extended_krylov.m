function [X,k,flag] = lr_lyap_extended_krylov(A,C1,itmax,tol)
%   [X,K,FLAG] = LR_LYAP_EXTENDED_KRYLOV(A,C1,ITMAX,TOL) risolve
%   l'equazione di Sylvester A*X + X*A' = C1*C1' per mezzo dell'algoritmo 
%   proiettivo basato sugli Spazi di Krylov Estesi.
%   INPUT:
%       - A,B,C1: matrici dei coefficienti
%       - ITMAX: numero massimo di iterazioni
%       - TOL: valore soglia per la condizione di arresto
%   OUTPUT:
%       - X: soluzione dell'equazione
%       - K: numero di iterazioni svolte
%       - FLAG: vale 0 o 1 e indica se le iterazioni sono state arrestate 
%         perche' l'errore ha iniziato a crescere 


%Operazioni preliminari
[~,s] = size(C1);
[LA,UA] = lu(A);
H = zeros(2*s);
itmin = 20;
errprec = 0;
err = 0;
flag = 0;
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
    m = k-1;
    if k>1
        errprec = err;
        err = sqrt(2)*norm(Atilde(2*m*s+1:2*m*s+s,2*(m-1)*s+1:2*m*s)*Y(2*(m-1)*s+1:2*m*s,:),'fro')/(norm(Y,'fro')*2*normA);
        if (err < tol || (errprec < err && k > itmin))
            if (errprec < err && k > itmin)
                disp('Warning: iterazioni arrestate per probabili errori di round-off');
                flag = 1;
            end
            X = V*Y*V';
            return
        end
    end
    
    %Aggiornamento dello spazio di Krylov
    %Partizionamento di V e W per ottenere Vtilde e Wtilde
    Vplus = V(:,end-2*s+1:end-s);
    Vmin = V(:,end-s+1:end);
    Vtilde = [A*Vplus,UA\(LA\Vmin)];
    
    %Ortogonalizzazione di Vtilde risepetto a V e di Wtilde rispetto a W
    %Poiche' le colonne di V e W sono ortonormali conviene calcolare come
    %sotto invece che con:
    %   Vtilde = Vtilde-V*V'*Vtilde
    %   Wtilde = Wtilde-W*W'*Wtilde
    for i = 1:k
        H = V(:,((i-1)*2*s+1):(i*2*s))'*Vtilde;
        Vtilde = Vtilde - V(:,((i-1)*2*s+1):(i*2*s))*H;
    end
    
    %Aggiornamento di V e W
    [Vtilde,~]=qr(Vtilde,0);
    V = [V,Vtilde];
end
disp('Warning: raggiunto il numero massimo di iterazioni');
X = V(:,1:end-2*s)*Y*V(:,1:end-2*s)';