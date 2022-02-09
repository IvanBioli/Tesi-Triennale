function [X,error1,error2] = low_rank_sylv_extended_krylov_old(A,B,C1,C2,itmax)
%AGGIUNGERE DESCRIZIONE

%Operazioni preliminari
[~,s]=size(C1);
[LA,UA]=lu(A);
[LB,UB]=lu(B');
error1 = zeros(itmax,1);
error2 = zeros(itmax,1);
%%%%%%%%%%%%%%%%%%%%% Corpo della procedura %%%%%%%%%%%%%%%%%%%%%%%%%
[V,~]=qr([A*C1,UA\(LA\C1)],0);
[W,~]=qr([B'*C2,UB\(LB\C2)],0);
for k = 1:itmax
    %Risoluzione del problema di taglia ridotta
    %Anche Atilde, e Btilde dovrebbero essere di Hessenberg, conviene
    %sfruttare questo nel calcolo?
    Atilde = V'*A*V;
    Btilde = W'*B*W;
    Ctilde = (V'*C1)*(W'*C2)';
    Y = sylv_bartels_stewart(Atilde,Btilde,Ctilde);
    %Controllo della convergenza, DA RIVEDERE
    %Xtilde = V*Y*W';
    %error1(k) = norm(A*Xtilde+Xtilde*B-C1*C2','fro')/norm(Y);
    m = k-1;
    if k>1
        error2(k) = sqrt(norm(Atilde(2*m*s+1:2*m*s+s,2*(m-1)*s+1:2*m*s)*Y(2*(m-1)*s+1:2*m*s,:),'fro')^2+norm(Btilde(2*(m-1)*s+1:2*m*s,2*m*s+1:2*m*s+s)'*Y(:,2*(m-1)*s+1:2*m*s)','fro')^2)/norm(Y);
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
        HV(((i-1)*2*s+1):(i*2*s),((k-1)*2*s+1):(k*2*s)) = V(:,((i-1)*2*s+1):(i*2*s))'*Vtilde;
        Vtilde = Vtilde - V(:,((i-1)*2*s+1):(i*2*s))*HV(((i-1)*2*s+1):(i*2*s),((k-1)*2*s+1):(k*2*s));
        HW(((i-1)*2*s+1):(i*2*s),((k-1)*2*s+1):(k*2*s)) = W(:,((i-1)*2*s+1):(i*2*s))'*Wtilde;
        Wtilde = Wtilde - W(:,((i-1)*2*s+1):(i*2*s))*HW(((i-1)*2*s+1):(i*2*s),((k-1)*2*s+1):(k*2*s));
    end
    %Aggiornamento di V e W
    [Vtilde,HV((k*2*s+1):((k+1)*2*s),((k-1)*2*s+1):(k*2*s))]=qr(Vtilde,0);
    V = [V,Vtilde];
    [Wtilde,HW((k*2*s+1):((k+1)*2*s),((k-1)*2*s+1):(k*2*s))]=qr(Wtilde,0);
    W = [W,Wtilde];
end
X =V*Y*W';