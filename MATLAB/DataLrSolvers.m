%Collezionamento di dati per gli algoritmi per solutori di rango basso per
%equazione di Sylvester e equazione di Lyapunov.

%Preliminari
n = 2000;           %Taglia della matrice A
m = 3000;            %Taglia della matrice B
s = 5;              %Colonne di C1 e C2
itmax = 100;        %Numero massimo di iterazioni
tol = 1e-13;         %Tolleranza per la condizione di arresto algoritmi non ADI
tolADI = 5e-13;         %Tolleranza per la condizione di arresto algoritmo ADI
N = 5;              %Numero di test
temp = 0;           %flag per decidere quali casi calcolare
NADI = 20;          %Numero parametri ADI

%Matrici per memorizzazione dati, una riga per ogni test. Nelle colonne:
%   1. Numero di iterazioni svolte
%   2. Tempo di esecuzione totale
%   3. Tempo di esecuzione di una iterazione
%   4. Errore relativo
%   5. Flag per indicare se le iterazioni sono state arrestate perchè
%      l'errore ha iniziato ad aumentare
%SylvData(:,:,i) è relativo a:
%   -i=1: lr_sylv_extended_krylov
%   -i=2: lr_sylv_krylov
%   -i=3: algoritmo di Bartels-Stewart
%   -i=4: funzione sylv di MATLAB
%se temp = 0 non vengono calcolati i casi i=2,2
%LyapData(:,:,i) è relativo a:
%   -i=1: lr_lyap_extended_krylov
%   -i=2: CF_ADI
%   -i=3: algoritmo di Bartels-Stewart
%   -i=4: funzione sylvester di MATLAB
%se temp = 0 non viene calcolato il caso i=3
SylvData = zeros(N,5,4);
LyapData = zeros(N,5,4);

%Matrici per memorizzazione dati medi. Nelle righe:
%   1. Numero medio di iterazioni svolte 
%   2. Tempo medio di esecuzione totale
%   3. Tempo medio di esecuzione di una iterazione
%   4. Errore relativo medio
%   5. Percentuale di iterazioni arrestate perchè l'errore ha iniziato ad 
%      aumentare
%SylvMean(:,i) e LyapMean(:,i) come prima
SylvMean = zeros(5,4); 
LyapMean = zeros(5,4);

%Equazione di Sylvester
for i = 1:N
    A = rand(n);
    A = -A'*A;          %A definita negativa
    B = rand(m);     
    B = -B'*B;          %B definita negativa
    C1 = rand(n,s);
    C2 = rand(m,s);
    C = C1*C2';
    
    tic
    [X,it,flag] = lr_sylv_extended_krylov(A,B,C1,C2,itmax,tol);
    time = toc;
    err = norm(A*X+X*B-C,'fro')/(norm(X,'fro')*(norm(A,'fro')+norm(B,'fro')));
    SylvData(i,1,1) = it;
    SylvData(i,2,1) = time;
    SylvData(i,3,1) = time/it;
    SylvData(i,4,1) = err;
    SylvData(i,5,1) = flag;
    
    if (temp == 1)
        tic
        [X,it] = lr_sylv_krylov(A,B,C1,C2,itmax,tol);
        time = toc;
        err = norm(A*X+X*B-C,'fro')/(norm(X,'fro')*(norm(A,'fro')+norm(B,'fro')));
        SylvData(i,1,2) = it;
        SylvData(i,2,2) = time;
        SylvData(i,3,2) = time/it;
        SylvData(i,4,2) = err;
        
        tic
        X = sylv_bartels_stewart(A,B,C);
        time = toc;
        err = norm(A*X+X*B-C,'fro')/(norm(X,'fro')*(norm(A,'fro')+norm(B,'fro')));
        SylvData(i,2,3) = time;
        SylvData(i,4,3) = err;
    end
    tic
    X = sylvester(A,B,C);
    time = toc;
    err = norm(A*X+X*B-C,'fro')/(norm(X,'fro')*(norm(A,'fro')+norm(B,'fro')));
    SylvData(i,2,4) = time;
    SylvData(i,4,4) = err;
end
for i = 1:4
    for j = 1:5
        SylvMean(j,i) = sum(SylvData(:,j,i))/N;
    end
end

%Equazione di Lyapunov
for i = 1:N
    A = rand(n);
    A = -A'*A;          %A definita negativa
    C1 = rand(n,s);
    C = C1*C1';
    
    tic
    [X,it,flag] = lr_lyap_extended_krylov(A,C1,itmax,tol);
    time = toc;
    err = norm(A*X+X*A'-C,'fro')/(norm(X,'fro')*2*norm(A,'fro'));
    LyapData(i,1,1) = it;
    LyapData(i,2,1) = time;
    LyapData(i,3,1) = time/it;
    LyapData(i,4,1) = err;
    LyapData(i,5,1) = flag;
    
    tic
    [Z,it] = CF_ADI(A,C1,NADI,itmax,tolADI);
    time = toc;
    Y = Z*Z';
    err = norm(A*Y+Y*A'+C1*C1','fro')/(norm(Y,'fro')*2*norm(A,'fro'));
    LyapData(i,1,2) = it;
    LyapData(i,2,2) = time;
    LyapData(i,3,2) = time/it;
    LyapData(i,4,2) = err;
    
    if (temp == 1)
        tic
        X = lyap_bartels_stewart(A,C);
        time = toc;
        err = norm(A*X+X*A'-C,'fro')/(norm(X,'fro')*2*norm(A,'fro'));
        LyapData(i,2,3) = time;
        LyapData(i,4,3) = err;
    end
    
    tic
    X = lyap(A,-C);
    time = toc;
    err = norm(A*X+X*A'-C,'fro')/(norm(X,'fro')*2*norm(A,'fro'));
    LyapData(i,2,4) = time;
    LyapData(i,4,4) = err;
end
for i = 1:4
    for j = 1:5
        LyapMean(j,i) = sum(LyapData(:,j,i))/N;
    end
end
