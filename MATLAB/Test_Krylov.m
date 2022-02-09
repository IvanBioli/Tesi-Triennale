%KRYLOV TEST

%Parametri
n = 1000;           %Taglia della matrice A
m = 1000;           %Taglia della matrice B
s = 1;              %Colonne di C1 e C2
itmax = 100;        %Numero massimo di iterazioni
tol = 5e-14;         %Tolleranza per la condizione di arresto
flag = 1;           %Flag per indicare se calcolare o meno l'errore effettivo

%Generazione di A, B, C1 e C2

A = rand(n);
A = -A'*A;          %A definita negativa
B = rand(m);
B = -B'*B;          %B definita negativa
C1 = rand(n,s);
C2 = rand(m,s);
C = C1*C2';

%TEST
%Funzione sylvester di MATLAB
tic
X = sylvester(A,B,C);
time = toc;
err = norm(A*X+X*B-C,'fro')/(norm(X,'fro')*(norm(A,'fro')+norm(B,'fro')));
fprintf('Funzione sylvester di MATLAB:\n\t Tempo impiegato: %f \n\t Errore relativo: %e\n',time,err);


%Funzione lr_sylv_krylov
tic
X = lr_sylv_krylov(A,B,C1,C2,itmax,tol);
time = toc;
err = norm(A*X+X*B-C,'fro')/(norm(X,'fro')*(norm(A,'fro')+norm(B,'fro')));
fprintf('Funzione lr_sylv_krylov:\n\t Tempo impiegato: %f \n\t Errore relativo: %e\n',time,err);
%Plot errore e numero di iterazioni
[~,it, error1,error2] = lr_sylv_krylov_error(A,B,C1,C2,itmax,tol,flag);
fprintf('\t Numero di iterazioni: %d \n',it);
if (flag ~= 0)
    fig1 = figure();
    semilogy(error1);
    hold on;
    semilogy(error2);
    legend('Errore reale','Errore "calcolato velocemente"');
    title(sprintf('Velocità di convergenza dell''algoritmo lr\\_sylv\\_krylov'));
    hold off
else
    fig1 = figure();
    semilogy(error1);
    legend('Errore "calcolato velocemente"');
    title(sprintf('Velocità di convergenza dell''algoritmo lr\\_sylv\\_krylov'));
end


%Funzione lr_sylv_extended_krylov
tic
[X,~,~]= lr_sylv_extended_krylov(A,B,C1,C2,itmax,tol);
time = toc;
err = norm(A*X+X*B-C,'fro')/(norm(X,'fro')*(norm(A,'fro')+norm(B,'fro')));
fprintf('Funzione lr_sylv_extended_krylov:\n\t Tempo impiegato: %f \n\t Errore relativo: %e\n',time,err);
%Plot errore e numero di iterazioni
[~,it, error1,error2] = lr_sylv_extended_krylov_error(A,B,C1,C2,itmax,tol,flag);
fprintf('\t Numero di iterazioni: %d \n',it);

if (flag ~= 0)
    fig2 = figure();
    semilogy(error1);
    hold on;
    semilogy(error2);
    legend('Errore reale','Errore "calcolato velocemente"');
    title(sprintf('Velocità di convergenza dell''algoritmo lr\\_sylv\\_extended\\_krylov'));
    hold off
else
    fig2 = figure();
    semilogy(error1);
    legend('Errore "calcolato velocemente"');
    title(sprintf('Velocità di convergenza dell''algoritmo lr\\_sylv_\\extended\\_krylov'));
end

%Funzione lyap di MATLAB
tic
X = lyap(A,-C1*C1');
time = toc;
err = norm(A*X+X*A'-C1*C1','fro')/(norm(X,'fro')*2*norm(A,'fro'));
fprintf('Funzione lyap di MATLAB:\n\t Tempo impiegato: %f \n\t Errore relativo: %e\n',time,err);

%Funzione lr_lyap_extended_krylov
tic
[X,~,~] = lr_lyap_extended_krylov(A,C1,itmax,tol);
time = toc;
err = norm(A*X+X*A'-C1*C1','fro')/(norm(X,'fro')*2*norm(A,'fro'));
fprintf('Funzione lr_lyap_extended_krylov:\n\t Tempo impiegato: %f \n\t Errore relativo: %e\n',time,err);
%Plot errore e numero di iterazioni
[~,it, error1,error2] = lr_lyap_extended_krylov_error(A,C1,itmax,tol,flag);
fprintf('\t Numero di iterazioni: %d \n',it);

if (flag ~= 0)
    fig2 = figure();
    semilogy(error1);
    hold on;
    semilogy(error2);
    legend('Errore reale','Errore "calcolato velocemente"');
    title(sprintf('Velocità di convergenza dell''algoritmo lr\\_lyap\\_extended\\_krylov'));
    hold off
else
    fig2 = figure();
    semilogy(error1);
    legend('Errore "calcolato velocemente"');
    title(sprintf('Velocità di convergenza dell''algoritmo lr\\_lyap_\\extended\\_krylov'));
end