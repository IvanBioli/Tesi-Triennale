function X = lyap_bartels_stewart(A,C)
% X = LYAP_BARTELS_STEWART(A,C) risolve l'equazione di Sylvester AX + XA* = C usando l'algoritmo di Bartels-Stewart basato sulla decomposizione di Schur complessa
% A,C: matrici dei coefficienti
% X: soluzione
[n,~] = size(C);
[U,R] = schur(A','complex');
C = U'*C*U;
X = zeros(n,n);
X(1,:) = C(1,:)/(conj(R(1,1))*eye(n) + R);
for k = 2:n
    v = C(k,:) - R(1:k-1,k)'*X(1:k-1,:);
    X(k,:) = v/(R + conj(R(k,k))*eye(n));
end
X = U*X*U';