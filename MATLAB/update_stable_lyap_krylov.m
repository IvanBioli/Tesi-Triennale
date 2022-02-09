function deltaX = update_stable_lyap_krylov(A,deltaA,deltaC,X,tol,UA,VA,UC,SigmaC)
%   DELTAX = UPDATE_STABLE_LYAP_ADI(A,DELTAA,DELTAC,X,TOL,UA,VA,UC,SIGMAC)
%   calcola una correzione DELTAX di rango basso tale che, data X soluzione
%   di AX+XA'=C, allora (X+DELTAX) e' soluzione dell'equazione di Sylvester
%       (A+DELTAA)(X+DELTAX)+(X+DELTAX)(A+DELTAA)'=(C+DELTAC)
%   usando come solutore per l'equazione di rango basso quello basato sul
%   metodo ADI.
%   INPUT:
%       -A: matrice dei coefficienti dell'equazione originaria
%       -DELTAA,DELTAC: matrici dei coefficienti della correzione
%       -X: soluzione dell'equazione originaria
%       -tol: tolleranza per la risoluzione dell'equazione della correzione
%       -UA,VA,UC,SIGMAC (OPZIONALI): matrici di fattorizzazione tali che
%        DELTAA = UA*VA', DELTAC = UC*SIGMAC*UC' 
%   OUTPUT:
%       -DELTAX: correzione

%Parametri e operazioni preliminari
toltsvd1 = 1e-7;
kA = -1;
kC = -1;
itmax = 1000;
N = 20;                     

%Se vengono date in input deltaA e deltaC ma non le loro fattorizzazioni, 
%allora UA,VA,UC,SigmaC vengono calcolate con la tsvd
if (nargin == 5)
    [UA,SigmaA,VA] = tsvd(deltaA,toltsvd1,kA);
    UA = UA*sqrt(SigmaA);
    VA = VA*sqrt(SigmaA);
    [UC,SigmaC,~] = tsvd(deltaC,toltsvd1,kC);
end

%Partizionamento del membro di destra della equazione per il calcolo della 
%correzione
[sC,~] = size(SigmaC);
[~,sA] = size(UA);
Utilde = [UC,UA,X*VA];
Sigma = [SigmaC,zeros(sC,2*sA);zeros(sA,sC+sA),-eye(sA);zeros(sA,sC),-eye(sA),zeros(sA,sA)];
[Qtilde,Rtilde] = qr(Utilde,0);
[Q,D] = eig(Rtilde*Sigma*Rtilde');
%Suddivisione in due equazioni di taglia ridotta
indexplus = find(diag(D)>0);
indexmin = find(diag(D)<0);
D1 = D(indexplus,indexplus);
Q1 = Q(:,indexplus);
D2 = -D(indexmin,indexmin);
Q2 = Q(:,indexmin);
U1 = Qtilde*Q1*sqrt(D1);
U2 = Qtilde*Q2*sqrt(D2);
A = A+deltaA;

%Calcolo della correzione
X1 = lr_lyap_extended_krylov(A,U1,itmax,tol);
X2 = lr_lyap_extended_krylov(A,U2,itmax,tol);
deltaX = X1-X2;