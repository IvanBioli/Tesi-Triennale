%{
v = [1,2,5,10,15];
s = 1;

n = 300;
DataCARE;
V1(1:N,:,1) = DataV1;
V2(1:N,:,1) = DataV2;
LR_ADI(1:N,:,1) = DataLR_ADI;
TotV1(1,:) = MeanV1;
TotV2(1,:) = MeanV2;
TotLR_ADI(1,:) = MeanLR_ADI;


n = 200;
DataCARE;
V1(1:N,:,2) = DataV1;
V2(1:N,:,2) = DataV2;
LR_ADI(1:N,:,2) = DataLR_ADI;
TotV1(2,:) = MeanV1;
TotV2(2,:) = MeanV2;
TotLR_ADI(2,:) = MeanLR_ADI;


n = 300;
DataCARE;
V1(1:N,:,3) = DataV1;
V2(1:N,:,3) = DataV2;
LR_ADI(1:N,:,3) = DataLR_ADI;
TotV1(3,:) = MeanV1;
TotV2(3,:) = MeanV2;
TotLR_ADI(3,:) = MeanLR_ADI;


n = 400;
DataCARE;
V1(1:N,:,4) = DataV1;
V2(1:N,:,4) = DataV2;
LR_ADI(1:N,:,4) = DataLR_ADI;
TotV1(4,:) = MeanV1;
TotV2(4,:) = MeanV2;
TotLR_ADI(4,:) = MeanLR_ADI;


n = 500;
DataCARE;
V1(1:N,:,5) = DataV1;
V2(1:N,:,5) = DataV2;
LR_ADI(1:N,:,5) = DataLR_ADI;
TotV1(5,:) = MeanV1;
TotV2(5,:) = MeanV2;
TotLR_ADI(5,:) = MeanLR_ADI;
%}

fig1 = figure();
plot(v,TotV1(:,1),'bo-')
hold on
plot(v,TotV2(:,1),'rs-')
plot(v,TotLR_ADI(:,1),'gd-')
legend('newton\_care\_v1','newton\_care\_v2','newton\_care\_ADI');
title('Tempo di esecuzione di un istanza');
hold off 

fig2 = figure();
semilogy(v,TotV1(:,2),'bo-')
hold on
semilogy(v,TotV2(:,2),'rs-')
semilogy(v,TotLR_ADI(:,2),'gd-')
legend('newton\_care\_v1','newton\_care\_v2','newton\_care\_ADI');
title('Plot in scala semilogaritmica dell''errore relativo commesso');
hold off 

fig3 = figure();
plot(v,TotV1(:,3),'bo-')
hold on
plot(v,TotV2(:,3),'rs-')
plot(v,TotLR_ADI(:,3),'gd-')
legend('newton\_care\_v1','newton\_care\_v2','newton\_care\_ADI');
title('Numero di iterazioni');
hold off 
