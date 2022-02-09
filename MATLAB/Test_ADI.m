%CF_ADI TEST

%Parametri
n = 1000;           %Taglia della matrice A
s = 1;              %Colonne di C1
N = 20;             %Numero dei parametri di shift calcolati per il metodo ADI
itmax = 250;        %Numero massimo di iterazioni
tol = 5e-13;         %Tolleranza per la condizione di arresto
flag = 1;           %Flag per indicare se calcolare o meno l'errore effettivo
density = 0.01;     %Densità (nel caso si usi sprandsym)
rc = 0.7;           %Reciproco del numero di condizionamento (nel caso si usi sprandsym)

%{
%Generazione di A e C1
A = rand(n);
A = -A'*A;                 %A definita negativa
C1 = rand(n,s);
C = C1*C1';
%}
%TEST
%Funzione lyap di MATLAB
tic
X = lyap(A,C);
time = toc;
err = norm(A*X+X*A'+C1*C1','fro')/(norm(X,'fro')*2*norm(A,'fro'));
fprintf('Funzione lyap di MATLAB:\n\t Tempo impiegato: %f \n\t Errore relativo: %e\n',time,err);

%Funzione CF_ADI
tic
Z = CF_ADI(A,C1,N,itmax,tol);
time = toc;
Y = Z*Z';
err = norm(A*Y+Y*A'+C1*C1','fro')/(norm(Y,'fro')*2*norm(A,'fro'));
fprintf('Funzione CF_ADI:\n\t Tempo impiegato: %f \n\t Errore relativo: %e\n',time,err);
%Plot errore e numero di iterazioni
[~,it, error1,error2] = CF_ADI_error(A,C1,N,itmax,tol,flag);
fprintf('\t Numero di iterazioni: %d \n',it);
if (flag ~= 0)
    semilogy(error1);
    hold on;
end
semilogy(error2);
if (flag ~= 0)
    legend('Errore reale','Errore "calcolato velocemente"');
else
    legend('Errore "calcolato velocemente"');
end