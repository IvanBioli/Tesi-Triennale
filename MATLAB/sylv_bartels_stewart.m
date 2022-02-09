function X = sylv_bartels_stewart(A,B,C)
% X=SYLV_BARTELS_STEWART(A,B,C) risolve l'equazione di Sylvester AX + XB = C usando l'algoritmo di Bartels-Stewart basato sulla decomposizione di Schur complessa
% A, B, C: matrici dei coefficienti
% X : soluzione
[n,m] = size(C);
[U,R] = schur(A','complex');
R = R';
[V,S] = schur(B,'complex');
C = U'*C*V;
X = zeros(n,m);
X(1,:) = C(1,:)/(R(1,1)*eye(m) + S);
for k = 2:n
    v = C(k,:) - R(k,1:k-1)*X(1:k-1,:);
    X(k,:) = v/(S + R(k,k)*eye(m));
end
X = U*X*V';