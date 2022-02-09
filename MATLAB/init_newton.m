function X = init_newton(A,B)
%   X=INIT_NEWTON(A,B) computes a stabilizing initial approximation X
%   for Newton’s method applied to the CARE C + XA + A’X - XBX = 0
%   A, B: matrix coefficients of the CARE
%   X: stabilizing initial approximation
n = size(A,1);
[U,TA] = schur(A);
TD = U'*B;
b = -min(real(ordeig(TA)));
b = max(0,b) + 0.5;
Z = lyapunov(TA+b*eye(n),2*TD*TD');
X = TD'*pinv(Z)*U';
if norm (X-X') > 1e-13
    M = X'*B*X+0.5*eye(n);
    X = lyapunov ((A-B*X)', -M);
end