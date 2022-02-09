function [Q,H] = ArnoldiInv(A,r,m)
%   [Q,H] = ARNOLDIINV(A,r,M) esegue M iterazioni del processo di Arnoldi a
%   partire dalla matrice A^(-1) N-by-N e da un vettore iniziale r (che non
%   necessariamente deve avere norma-2 unitaria). Per M < N produce in 
%   output un matrice Q N-by-(M+1) con colonne ortonormali e una matrice di
%   Hessenberg superiore H (M+1)-by-M tale che:
%       A^(-1)*Q(:,1:M) = Q(:,1:M)*H(1:M,1:M) + H(M+1,M)*Q(:,M+1)*e_M'
%   dove e_M e' l'M-esima colonna della matrice identita' di taglia M.

[n,~]=size(A);
r = r/norm(r);
Q = zeros(n,m); 
Q(:,1) = r;
H = zeros(m+1,m);
tol = 1e-16;

for k=1:m
    z = A\Q(:,k);               %Non viene calcolata inv(A) ma si risolvono sitemi lineari
    for i=1:k
        H(i,k) = Q(:,i)'*z;
        z = z - H(i,k)*Q(:,i);
    end
    if k < n
       H(k+1,k) = norm(z);
       if (H(k+1,k) < tol)
           return
       end
       Q(:,k+1) = z/H(k+1,k);
   end
end