function [X,error1,error2] = low_rank_sylv_krylov_old(A,B,C1,C2,itmax)
%AGGIUNGERE DESCRIZIONE

%%%%%%%%%%%%%%%%%%%%% Corpo della procedura %%%%%%%%%%%%%%%%%%%%%%%%%
[~,s]=size(C1);
[V,~]=qr(A*C1,0);
[W,~]=qr(B'*C2,0);
C = C1*C2';
for k = 1:itmax
    %Risoluzione del problema di taglia ridotta
    %Anche Atilde, e Btilde dovrebbero essere di Hessenberg, conviene
    %sfruttare questo nel calcolo?
    Atilde = V'*A*V;
    Btilde = W'*B*W;
    C1tilde = V'*C1;
    C2tilde = W'*C2;
    Ctilde = C1tilde*C2tilde';
    if (k > 1)
        error2(k-1) = norm(HV*Y*[eye(k-1),zeros(k-1,s)]+[eye(k-1);zeros(s,k-1)]*Y*HW'-Ctilde,'fro')/norm(Y,'fro');
    end
    Y = sylv_bartels_stewart(Atilde,Btilde,Ctilde);
    %Controllo della convergenza, DA RIVEDERE
    X = V*Y*W';
    error1(k) = norm(A*X+X*B-C,'fro')/norm(Y,'fro');
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