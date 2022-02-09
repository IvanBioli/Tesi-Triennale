function p = ADI_Suboptimal(A,c,kplus,kmin)
%   P = ADI_SUBOPTIMAL(A,C,KPLUS,KMIN) calcola dei parametri subottimali
%   per il metodo ADI secondo l'algoritmo euristico di T. Penzl.
%   Applicabile a matrici definite negative oppure tali che R sia contenuto
%   nel semipiano complesso sinistro.
%   INPUT:
%       - A: matrice dei coefficienti dell'equazione di Lyapunov
%       - C: numero di parametri da ottenere
%       - KPLUS: numero di iterazioni del metodo di Arnoldi per A 
%       - KMIN: numero di iterazioni del metodo di Arnoldi per A^(-1)
%   Se KPLUS+KMIN<C vengono impostati di defalut KPLUS = KPLUS+C e 
%   KMIN = KMIN +C. 
%   OUTPUT:
%       - P: vettore contenente i parametri di shift

%Operazioni preliminari
[n,~]=size(A);
count = 0;
k = kplus+kmin;
p = zeros(c,1);
if (k<c)
    kplus = kplus+c;
    kmin = kmin+c;
    k = k+2*c;
    disp('Warning: kplus+kmin < c, impostato di defalut kplus = kplus+c e kmin = kmin+c');
end

%Definizione della funzione ausiliaria s
s = @(v,t) abs(prod(v-t)/prod(v+t));

%Calcolo dell'insieme R
r = rand(n,1);
[~,H] = Arnoldi(A,r,kplus);
[~,W] = ArnoldiInv(A,r,kmin);
Rplus = eig(H(1:kplus,1:kplus));
Rmin = 1./eig(W(1:kmin,1:kmin));
R = [Rplus;Rmin];
%Verifica che R sia contenuto nel semipiano sinistro 
if ~all(real(R)<0)
    error('L''insieme R non e'' contenuto nel semipiano complesso sinistro');
end

%Passo iniziale
aux1 = zeros(k,k);
for i = 1:k
    for j = 1:k
        aux1(i,j) = abs((R(i)-R(j))/(R(i)+R(j)));
    end
end
[~,i] = min(max(aux1));
count = count+1;
p(count) = R(i);
R(R == p(count)) = [];
if ~isreal(p(count))
    count = count+1;
    p(count) = conj(p(count-1));
    R(R == p(count)) = [];
end

%Passi successivi al primo
while (count < c)
    imax = 1;
    val = s(p,R(1));
    for i = 2:k-count
        temp = s(p,R(i));
        if (temp > val)
            imax = i;
            val = temp;
        end
    end
    count = count+1;
    p(count) = R(imax);
    R(R == p(count)) = [];
    if ~isreal(p(count))
        count = count+1;
        p(count) = conj(p(count-1));
        R(R == p(count)) = [];
    end
end