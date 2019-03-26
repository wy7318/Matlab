function [Q,R,p] = mymodifiedgramschmidtqr(A)
% Modified Gram-Schmidt QR with column pivoting.
% Q is orthogonal, U is upper triangular
% The relation A(:,p) = Q*R holds (up to round-off errors)
[m,n] = size(A);
R = zeros(n,n);
p = 1:n;
for k = 1:n
    v = zeros(1,n); for j = k:n, v(j) = norm(A(:,j)); end
    [vmax,jmax] = max(v(k:n)); % vmax is not used
    jm = jmax(1); % there may be more than one value
    k1 = k+[0,jm-1]; k2 = k+[jm-1,0];
    A(:,k1) = A(:,k2); p(k1) = p(k2); % permute
    R(k,k) = norm(A(:,k));
    A(:,k) = A(:,k)/R(k,k);
    i = k+1:n;
    R(k,i) = A(:,k)'*A(:,i);
    A(:,i) = A(:,i)-A(:,k)*R(k,i);
end
Q = A;
