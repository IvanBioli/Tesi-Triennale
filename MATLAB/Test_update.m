%TEST PER GLI ALGORITMI DI UPDATE
%Parametri
n = 1000;           %Taglia della matrice A
m = 1000;           %Taglia della matrice B
sA = 1;             %Rango di deltaA
sB = 1;             %Rango di deltaB
sC = 5;             %Rango di deltaC
itmax = 1000;       %Numero massimo di iterazioni
tol = 1e-7;         %Tolleranza per la condizione di arresto
flag = 1;           %Flag per indicare se calcolare o meno l'errore effettivo
density = 0.01;     %Densità (nel caso si usi sprandsym)
rc = 0.7;           %Reciproco del numero di condizionamento (nel caso si usi sprandsym)


%TEST PER L'ALGORITMO DI UPDATE PER L'EQUAZIONE DI SYLVESTER BASATO SUI
%SOTTOSPAZI DI KRYLOV
%Matrici dei coefficienti
A = rand(n);
A = -A'*A;          %A definita negativa
UA = rand(n,sA);
VA = rand(n,sA);
deltaA = UA*VA';
B = rand(m);     
B = -B'*B;          %B definita negativa
UB = rand(m,sB);
VB = rand(m,sB);
deltaB = UB*VB';
C = rand(n,m);
UC = rand(n,sC);
VC = rand(m,sC);
deltaC = UC*VC';

%Caclolo della soluzione e dell'update
tic
X = sylv_bartels_stewart(A,B,C);
time = toc;
err1 = norm(A*X+X*B-C,'fro')/norm(X,'fro');
fprintf('Equazione di Sylvester, algoritmo di Bartels-Stewart:\n\t Tempo impiegato: %e \n\t Errore: %e\n',time,err1);
tic
deltaX = update_sylv_krylov(A,deltaA,B,deltaB,deltaC,X,tol,UA,VA,UB,VB,UC,VC);
time = toc;
X = X+deltaX;
err2 = norm((A+deltaA)*X+X*(B+deltaB)-(C+deltaC),'fro')/norm(X,'fro');
fprintf('Equazione di Sylvester, solutore basato sui sottospazi di Krylov Estesi:\n\t Tempo impiegato: %e \n\t Errore originario: %e \n\t Errore soluzione updated: %e\n',time,err1,err2);

tic
deltaX = update_sylv_krylov(A+deltaA,deltaA,B+deltaB,deltaB,deltaC,X,tol,UA,VA,UB,VB,UC,VC);
time = toc;
X = X+deltaX;
err3 = norm((A+2*deltaA)*X+X*(B+2*deltaB)-(C+2*deltaC),'fro')/norm(X,'fro');
fprintf('Equazione di Sylvester, solutore basato sui sottospazi di Krylov Estesi, secondo update:\n\t Tempo impiegato: %e \n\t Errore originario: %e \n\t Errore soluzione updated: %e\n',time,err2,err3);

%TEST PER L'ALGORITMO DI UPDATE PER L'EQUAZIONE DI LYAPUNOV STABILE
%Matrici dei coefficienti
A = rand(n);
A = -A'*A;          %A definita negativa
UA = rand(n,sA);
%VA = rand(n,sA);
VA = -UA;
deltaA = UA*VA';
C = rand(n);
C = -C'*C;          %C definita negativa
UC = rand(n,sC);
SigmaC = -diag(rand(sC,1));
deltaC = UC*SigmaC*UC';

%Caclolo della soluzione e dell'update
tic
X = lyap_bartels_stewart(A,C);
time = toc;
err1 = norm(A*X+X*A'-C,'fro')/norm(X,'fro');
fprintf('Equazione di Lyapunov, algoritmo di Bartels-Stewart:\n\t Tempo impiegato: %e \n\t Errore: %e\n',time,err1);
tic
deltaX1 = update_stable_lyap_ADI(A,deltaA,deltaC,X,tol,UA,VA,UC,SigmaC);
time1 = toc;
tic
deltaX2 = update_stable_lyap_krylov(A,deltaA,deltaC,X,tol,UA,VA,UC,SigmaC);
time2 = toc;
X1 = X+deltaX1;
A = A+deltaA;
err2 = norm(A*X1+X1*A'-(C+deltaC),'fro')/norm(X1,'fro');
fprintf('Equazione di Lyapunov, solutore basato sul metodo CF_ADI:\n\t Tempo impiegato: %e \n\t Errore originario: %e \n\t Errore soluzione updated: %e\n',time1,err1,err2);

X2 = X+deltaX2;
err3 = norm(A*X2+X2*A'-(C+deltaC),'fro')/norm(X2,'fro');
fprintf('Equazione di Lyapunov, solutore basato sul metodo di Krylov:\n\t Tempo impiegato: %e \n\t Errore originario: %e \n\t Errore soluzione updated: %e\n',time2,err1,err3);
