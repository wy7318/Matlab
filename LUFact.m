function [L,U] = mylu(A)
n = size(A,1);
for k = 1:n
    if A(k,k)==0
        warning('LU factorization fails');
        L = []; U = []; return; 
    end
    i = k+1:n;
    A(i,k) = A(i,k)/A(k,k);
    A(i,i) = A(i,i)-A(i,k)*A(k,i);
end
L = tril(A,-1)+eye(n); U = triu(A); 
