%Collezionamento di dati per gli algoritmi per il metodo di Newton
%n = 100;                    %Taglia della matrice
%s = 1;                      %Colonne di BU
tol = 5e-13;                 %Tolleranza
itmax = 40;                %Numero massimo di iterazioni
N = 5;                      %Numero di test da effettuare

DataV1 = zeros(N,3);
DataV2 = zeros(N,3);
DataLR_krylov = zeros(N,3);
DataLR_ADI = zeros(N,3);

for i = 1:N
    disp(i)
    A = rand(n);
    A = -A'*A;
    BU = rand(n,s);
    B = BU*(BU)';
    C = rand(n);
    C = -C*C';
    X0 = zeros(n);

    tic
    [X,it] = newton_care_v1(A,B,C,X0,tol,itmax);
    DataV1(i,1) = toc;
    DataV1(i,2) = norm(X*A + A'*X - X*B*X - C,'fro')/(norm(X,'fro')*(norm(A,'fro')+norm(B,'fro')));
    DataV1(i,3) = it;
    
    tic
    [X,it] = newton_care_v2(A,B,C,X0,tol,itmax);
    DataV2(i,1) = toc;
    DataV2(i,2) = norm(X*A + A'*X - X*B*X - C,'fro')/(norm(X,'fro')*(norm(A,'fro')+norm(B,'fro')));
    DataV2(i,3) = it;
    %{
    try
        tic
        [X,it] = lr_newton_care_krylov(A,BU,C,X0,tol,itmax);
        DataLR_krylov(i,1) = toc;
        DataLR_krylov(i,2) = norm(X*A + A'*X - X*B*X - C);
        DataLR_krylov(i,3) = it;
    catch
        DataLR_krylov(i,1) = NaN;
        DataLR_krylov(i,2) = NaN;
        DataLR_krylov(i,2) = NaN;
    end
    %}
    try
        tic
        [X,it] = lr_newton_care_ADI(A,BU,C,X0,tol,itmax);
        DataLR_ADI(i,1) = toc;
        DataLR_ADI(i,2) = norm(X*A + A'*X - X*B*X - C,'fro')/(norm(X,'fro')*(norm(A,'fro')+norm(B,'fro')));
        DataLR_ADI(i,3) = it;
    catch
        DataLR_ADI(i,1) = NaN;
        DataLR_ADI(i,2) = NaN;
        DataLR_ADI(i,2) = NaN;
    end
end

MeanV1 = nanmean(DataV1);
MeanV1(4) = 0;
MeanV2 = nanmean(DataV2);
MeanV2(4) = 0; 
MeanLR_krylov = nanmean(DataLR_krylov);
nan = sum(isnan(DataLR_krylov));
MeanLR_krylov(4) = nan(1)/N; 
MeanLR_ADI = nanmean(DataLR_ADI);
nan = sum(isnan(DataLR_ADI));
MeanLR_ADI(4) = nan(1)/N; 