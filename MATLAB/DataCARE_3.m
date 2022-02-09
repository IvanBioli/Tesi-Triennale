%{
n = 400;
temp = n/100;
k = 2;
DataCARE;
V1(k,1:3,temp) = DataV1;
V2(k,1:3,temp) = DataV2;
LR_ADI(k,1:3,temp) = DataLR_ADI;
LR_ADI(k,4,temp) = 1;
%}
for temp = 1:4
TotV1(temp,1:3) = nanmean(V1(1:5,1:3,temp));
TotV2(temp,1:3) = nanmean(V2(1:5,1:3,temp));
TotLR_ADI(temp,1:3) = nanmean(LR_ADI(1:5,1:3,temp));
TotLR_ADI(temp,4) = sum(isnan(LR_ADI(1:5,3,temp)))/5;
end