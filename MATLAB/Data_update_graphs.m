%Grafici aggiuntivi per update
sylvbart = zeros(upd+2,1);
sylvkryl = zeros(upd+2,1);
lyapbart = zeros(upd+2,1);
lyapkryl = zeros(upd+2,1);
lyapadi = zeros(upd+2,1);

for i=1:(upd+1)
    sylvbart(i+1) = sum(DataSylvBartels(1,1:i,1));
    sylvkryl(i+1) = sum(DataSylvKrylov(1,1:i,1));
    lyapbart(i+1) = sum(DataLyapBartels(1,1:i,1));
    lyapadi(i+1) = sum(DataLyapADI(1,1:i,1));
    lyapkryl(i+1) = sum(DataLyapKrylov(1,1:i,1));
end

fig5 = figure();
plot(0:(upd+1),sylvbart,'bo-');
hold on
plot(0:(upd+1),sylvkryl,'rs-');
legend('Algoritmo di Bartels-Stewart','Algoritmo di update di Krylov');
title('Plot del tempo di esecuzione totale dopo ogni passo di update');
hold off 

fig6 = figure();
plot(0:(upd+1),lyapbart,'bo-');
hold on
plot(0:(upd+1),lyapkryl,'rs-');
plot(0:(upd+1),lyapadi,'gd-');
legend('Algoritmo di Bartels-Stewart','Algoritmo di update di Krylov','Algoritmo di update ADI');
title('Plot del tempo di esecuzione totale dopo ogni passo di update');
hold off 