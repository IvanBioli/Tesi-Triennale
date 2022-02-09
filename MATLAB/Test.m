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

