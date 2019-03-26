function [Q,R,p] = mygivensqr(A)
% Givens QR with column pivoting.
% Q is orthogonal, U is upper triangular
% The relation A(:,p) = Q*R holds (up to round-off errors)
[m,n] = size(A);
Q = speye(m,m);
p = 1:n; % initialize permutation vector
for k = 1:min(m-1,n)
    i = k:m;
    v = zeros(1,n); for j = k:n, v(j) = norm(A(i,j)); end
    [vmax,jmax] = max(v(k:n)); % vmax is not used
    jm = jmax(1); % there may be more than one value
    k1 = k+[0,jm-1]; k2 = k+[jm-1,0];
    A(:,k1) = A(:,k2); p(k1) = p(k2); % permute
     j = k:n;
     for i = k+1:m
        G = planerot(A([k,i],k)); % Givens matrix
        A([k,i],j) = G*A([k,i],j); % update A
        Q([k,i],:) = G*Q([k,i],:); % update Q
    end
end
R = triu(A); % extract R
k = 1:min(m,n); Q = full(Q(k,:)'); R = R(k,:); % lean Q and R
