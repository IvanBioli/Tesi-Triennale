function deltaX = update_sylv_krylov(A,deltaA,B,deltaB,deltaC,X,tol,UA,VA,UB,VB,UC,VC)
%   DELTAX = UPDATE_SYLV_KRYLOV(A,DELTAA,B,DELTAB,DELTAC,X,TOL,UA,VA,UB,VB,UC,VC)
%   calcola una correzione DELTAX di rango basso tale che, data X soluzione
%   di AX+XB=C, allora (X+DELTAX) e' soluzione dell'equazione di Sylvester
%       (A+DELTAA)(X+DELTAX)+(X+DELTAX)(B+DELTAB)=(C+DELTAC)
%   usando come solutore per l'equazione di rango basso quello basato sugli
%   Spazi di Krylov Estesi.
%   INPUT:
%       -A,B: matrici dei coefficienti dell'equazione originaria
%       -DELTAA,DELTAB,DELTAC: matrici dei coefficienti della correzione
%       -X: soluzione dell'equazione originaria
%       -tol: tolleranza per la risoluzione dell'equazione della correzione
%       -UA,VA,UB,VB,UC,VC (OPZIONALI): matrici di fattorizzazione tali che
%        DELTAA = UA*VA', DELTAB = UB*VB' e DELTAC = UC*VC' 
%   OUTPUT:
%       -DELTAX: correzione

%Parametri e operazioni preliminari
toltsvd1 = 1e-7;
toltsvd2 = 1e-7;
kA = -1;
kB = -1;
kC = -1;
kUV = -1;
itmax = 500;

%Se vengono date in input deltaA,deltaB e deltaC ma non le loro 
%fattorizzazioni, allora UA,VA,UB,VB,UC,VC vengono calcolate con la tsvd
if (nargin == 7)
    [UA,SigmaA,VA] = tsvd(deltaA,toltsvd1,kA);
    UA = UA*sqrt(SigmaA);
    VA = VA*sqrt(SigmaA);
    [UB,SigmaB,VB] = tsvd(deltaB,toltsvd1,kB);
    UB = UB*sqrt(SigmaB);
    VB = VB*sqrt(SigmaB);
    [UC,SigmaC,VC] = tsvd(deltaC,toltsvd1,kC);
    UC = UC*sqrt(SigmaC);
    VC = VC*sqrt(SigmaC);
end

%Generazione e successiva compressione di U e V tali che:
%   deltaC-deltaA*X-X*deltaB = U*V'
U = [UC,-UA,-X*UB];
V = [VC,X'*VA,VB];
[QU,RU] = qr(U,0);
[QV,RV] = qr(V,0);
[UR,SigmaR,VR] = tsvd(RU*RV',toltsvd2,kUV);
U = QU*UR*sqrt(SigmaR);
V = QV*VR*sqrt(SigmaR);

%Risoluzione dell'equazione per la correzione con il metodo proiettivo
%basato su Spazi di Krylov Estesi
A = A+deltaA;
B = B+deltaB;
deltaX = lr_sylv_extended_krylov(A,B,U,V,itmax,tol);