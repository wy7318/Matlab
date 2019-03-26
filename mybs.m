function x = mybs(U,y)
n = size(U,2);
x = zeros(n,1);
for i = n:-1:1
    if U(i,i)==0
        warning('U has 0 diagonal coef');
        x = []; return
    end
    x(i) = (y(i)-U(i,i+1:n)*x(i+1:n))/U(i,i);
end