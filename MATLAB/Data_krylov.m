%Collezionamento di dati per gli algoritmi lr_sylv_extended_krylov e 
%lr_lyap_extended_krylov 

%Preliminari
n = 1000;           %Taglia della matrice A
m = 1000;           %Taglia della matrice B
s = 1;              %Colonne di C1 e C2
itmax = 1000;       %Numero massimo di iterazioni
tol = 1e-9;         %Tolleranza per la condizione di arresto
N = 1;             %Numero di test

%Matrici per memorizzazione dati, una riga per ogni test. Nelle colonne:
%   1. Numero di iterazioni svolte
%   2. Tempo di esecuzione totale
%   3. Tempo di esecuzione di una iterazione
%   4. Errore relativo
%   5. Flag per indicare se le iterazioni sono state arrestate perchè
%      l'errore ha iniziato ad aumentare
SylvData = zeros(N,5);
LyapData = zeros(N,5);
%Matrici per memorizzazione dati medi. Nelle colonne:
%   1. Numero medio di iterazioni svolte 
%   2. Tempo medio di esecuzione totale
%   3. Tempo medio di esecuzione di una iterazione
%   4. Errore relativo medio
%   5. Percentuale di iterazioni arrestate perchè l'errore ha iniziato ad 
%      aumentare
SylvMean = zeros(1,5); 
LyapMean = zeros(1,5);

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
    SylvData(i,1) = it;
    SylvData(i,2) = time;
    SylvData(i,3) = time/it;
    SylvData(i,4) = err;
    SylvData(i,5) = flag; 
end
for j = 1:5
    SylvMean(j) = sum(SylvData(:,j))/N;
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
    LyapData(i,1) = it;
    LyapData(i,2) = time;
    LyapData(i,3) = time/it;
    LyapData(i,4) = err;
    LyapData(i,5) = flag; 
end
for j = 1:5
    LyapMean(j) = sum(LyapData(:,j))/N;
end