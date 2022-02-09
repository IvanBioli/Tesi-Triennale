%Collezionamento di dati per l'algoritmo CF_ADI

%Preliminari
n = 1000;           %Taglia della matrice A
m = 1000;           %Taglia della matrice B
s = 1;              %Colonne di C1
itmax = 250;       %Numero massimo di iterazioni
tol = 5e-8;         %Tolleranza per la condizione di arresto
N = 1;              %Numero di test
NADI = 30;          %Numero parametri ADI

%Matrici per memorizzazione dati, una riga per ogni test. Nelle colonne:
%   1. Numero di iterazioni svolte
%   2. Tempo di esecuzione totale
%   3. Tempo di esecuzione di una iterazione
%   4. Errore relativo
%   5. Flag per indicare se le iterazioni sono state arrestate perchè
%      l'errore ha iniziato ad aumentare
Data = zeros(N,5); 
%Matrici per memorizzazione dati medi. Nelle colonne:
%   1. Numero medio di iterazioni svolte 
%   2. Tempo medio di esecuzione totale
%   3. Tempo medio di esecuzione di una iterazione
%   4. Errore relativo medio
%   5. Percentuale di iterazioni arrestate perchè l'errore ha iniziato ad 
%      aumentare
Mean = zeros(1,5); 


%Matrici non sparse
for i = 1:N
    A = rand(n);
    A = -A'*A;          %A definita negativa
    C1 = rand(n,s);
    C = C1*C1';
    
    tic
    [Z,it] = CF_ADI(A,C1,NADI,itmax,tol);
    time = toc;
    Y = Z*Z';
    err = norm(A*Y+Y*A'+C1*C1','fro')/norm(Y,'fro');
    Data(i,1) = it;
    Data(i,2) = time;
    Data(i,3) = time/it;
    Data(i,4) = err;
end
for j = 1:5
    Mean(j) = sum(Data(:,j))/N;
end