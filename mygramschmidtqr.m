function [Q,R,p] = mygramschmidtqr(A)
% Gram-Schmidt QR with column pivoting.
% Q is orthogonal, U is upper triangular
% The relation A(:,p) = Q*R holds (up to round-off errors)
[m,n] = size(A);
R = zeros(n,n);
v = zeros(1,n);
for j = 1:n, v(j) = norm(A(:,j)); end
[v,p] = sort(v,'descend');
A = A(:,p);
for j = 1:n
    k = 1:j-1;
    R(k,j) = A(:,k)'*A(:,j);
    A(:,j) = A(:,j)-A(:,k)*R(k,j);
    R(j,j) = norm(A(:,j));
    A(:,j) = A(:,j)/R(j,j);
end
Q = A;
