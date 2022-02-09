%Collezionamento di dati per gli algoritmi di update

%Parametri
n = 1000;           %Taglia della matrice A
m = 1000;           %Taglia della matrice B
sA = 1;             %Rango di deltaA
sB = 1;             %Rango di deltaB
sC = 1;             %Rango di deltaC
itmax = 1000;       %Numero massimo di iterazioni
tol = 1e-13;        %Tolleranza per la condizione di arresto di algoritmi non ADI
tolADI = 5e-13;     %Tolleranza per la condizione di arresto metodo ADI
flag = 1;           %Flag per indicare se calcolare o meno l'errore effettivo
density = 0.01;     %Densità (nel caso si usi sprandsym)
rc = 0.7;           %Reciproco del numero di condizionamento (nel caso si usi sprandsym)
N = 1;              %Numero di test da eseguire
upd = 5;            %Numero di update per ogni test


%EQUAZIONE DI SYLVESTER  
DataSylvBartels = zeros(N,upd+1,2);
DataSylvKrylov = zeros(N,upd+1,2);
%Esecuzione dei test
for i = 1:N
    %Risoluzione equazione iniziale
    A = rand(n);
    A = -A'*A;          %A definita negativa
    B = rand(m);     
    B = -B'*B;          %B definita negativa
    C = rand(n,m);
    tic
    X = sylv_bartels_stewart(A,B,C);
    DataSylvBartels(i,1,1) = toc;
    DataSylvKrylov(i,1,1) = DataSylvBartels(i,1,1);
    DataSylvBartels(i,1,2) = norm(A*X+X*B-C,'fro')/(norm(X,'fro')*(norm(A,'fro')+norm(B,'fro')));
    DataSylvKrylov(i,1,2) = DataSylvBartels(i,1,2);
    XKryl = X;
    XBart = X;
    %Updates
    for j = 2:upd+1
        UA = rand(n,sA);
        %VA = rand(n,sA);
        VA = -UA;
        deltaA = UA*VA';
        UB = rand(m,sB);
        %VB = rand(m,sB);
        VB = -UB;
        deltaB = UB*VB';
        UC = rand(n,sC);
        VC = rand(m,sC);
        deltaC = UC*VC';
        
        tic
        XBart = sylv_bartels_stewart(A+deltaA,B+deltaB,C+deltaC);
        DataSylvBartels(i,j,1) = toc;
        
        tic
        deltaXKryl = update_sylv_krylov(A,deltaA,B,deltaB,deltaC,XKryl,tol,UA,VA,UB,VB,UC,VC);
        XKryl = XKryl + deltaXKryl;
        DataSylvKrylov(i,j,1) = toc;
        
        A = A+deltaA;
        B = B+deltaB;
        C = C+deltaC;
        DataSylvBartels(i,j,2) = norm(A*XBart+XBart*B-C,'fro')/(norm(XBart,'fro')*(norm(A,'fro')+norm(B,'fro')));
        DataSylvKrylov(i,j,2) = norm(A*XKryl+XKryl*B-C,'fro')/(norm(XKryl,'fro')*(norm(A,'fro')+norm(B,'fro')));
    end
end

%Grafici
fig1 = figure();
semilogy(DataSylvBartels(1,:,2),'bo-');
hold on
semilogy(DataSylvKrylov(1,:,2),'rs-');
legend('Algoritmo di Bartels-Stewart','Algoritmo di update di Krylov');
title('Plot in scala semilogaritmica: errore per ogni passo di update');
hold off
%saveas(fig1,'sylv_error','png');

fig2 = figure();
plot(DataSylvBartels(1,:,1),'bo-');
hold on
plot(DataSylvKrylov(1,:,1),'rs-');
legend('Algoritmo di Bartels-Stewart','Algoritmo di update di Krylov');
title('Plot del tempo di esecuzione per ogni passo di update');
hold off 
%saveas(fig2,'sylv_time','png');


%EQUAZIONE DI LYAPUNOV
DataLyapBartels = zeros(N,upd+1,2);
DataLyapKrylov = zeros(N,upd+1,2);
DataLyapADI = zeros(N,upd+1,2);
%Esecuzione dei test
for i = 1:N
    %Risoluzione equazione iniziale
    A = rand(n);
    A = -A'*A;          %A definita negativa
    C = rand(n);
    C = -C'*C;          %C definita negativa
    
    tic
    Y = lyap_bartels_stewart(A,C);
    DataLyapBartels(i,1,1) = toc;
    DataLyapKrylov(i,1,1) = DataLyapBartels(i,1,1);
    DataLyapADI(i,1,1) = DataLyapBartels(i,1,1);
    DataLyapBartels(i,1,2) = norm(A*Y+Y*A'-C,'fro')/(norm(Y,'fro')*2*norm(A,'fro'));
    DataLyapKrylov(i,1,2) = DataLyapBartels(i,1,2);
    DataLyapADI(i,1,2) = DataLyapBartels(i,1,2);
    YKryl = Y;
    YBart = Y;
    YADI = Y;
    %Updates
    for j = 2:upd+1
        UA = rand(n,sA);
        %VA = rand(n,sA);
        VA = -UA;
        deltaA = UA*VA';
        UC = rand(n,sC);
        SigmaC = -diag(rand(sC,1));
        deltaC = UC*SigmaC*UC';
        
        tic
        YBart = lyap_bartels_stewart(A+deltaA,C+deltaC);
        DataLyapBartels(i,j,1) = toc;
        
        tic
        deltaYKryl = update_stable_lyap_krylov(A,deltaA,deltaC,YKryl,tol,UA,VA,UC,SigmaC);
        YKryl = YKryl + deltaYKryl;
        DataLyapKrylov(i,j,1) = toc;
        
        tic
        deltaYADI = update_stable_lyap_ADI(A,deltaA,deltaC,YADI,tolADI,UA,VA,UC,SigmaC);
        YADI = YADI + deltaYADI;
        DataLyapADI(i,j,1) = toc;

        A = A+deltaA;
        B = B+deltaB;
        C = C+deltaC;
        DataLyapBartels(i,j,2) = norm(A*YBart+YBart*A'-C,'fro')/(norm(YBart,'fro')*2*norm(A,'fro'));
        DataLyapKrylov(i,j,2) = norm(A*YKryl+YKryl*A'-C,'fro')/(norm(YKryl,'fro')*2*norm(A,'fro'));
        DataLyapADI(i,j,2) = norm(A*YADI+YADI*A'-C,'fro')/(norm(YADI,'fro')*2*norm(A,'fro'));
    end
end

%Grafici
fig3 = figure();
semilogy(DataLyapBartels(1,:,2),'bo-');
hold on
semilogy(DataLyapKrylov(1,:,2),'rs-');
semilogy(DataLyapADI(1,:,2),'gd-');
legend('Algoritmo di Bartels-Stewart','Algoritmo di update di Krylov','Algortimo di update ADI');
title('Plot in scala semilogaritmica: errore per ogni passo di update');
hold off
%saveas(fig3,'lyap_error','png');

fig4 = figure();
plot(DataLyapBartels(1,:,1),'bo-');
hold on
plot(DataLyapKrylov(1,:,1),'rs-');
plot(DataLyapADI(1,:,1),'gd-');
legend('Algoritmo di Bartels-Stewart','Algoritmo di update di Krylov','Algoritmo di update ADI');
title('Plot del tempo di esecuzione per ogni passo di update');
hold off 
%saveas(fig4,'lyap_time','png');