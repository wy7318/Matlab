function [x,cv] = mybasic(A,b,M,x,TOL,nmax)
% Basic iteration x <- x+M\(b-A*x)
x = x(:); % make sure initial guess is column vector
r = b-A*x;
e = M\r;
err0 = norm(e);  % euclidean magnitude
rerr = 1;        % relative error
n = 0;
if nargout==2
    cv = zeros(nmax+1,2+length(x)); % preallocate space for cv
    cv(1,:) = [n,x.',rerr]; % iter. #, iterate (row vector), rel. error
end
while (rerr>TOL && n<nmax)
    x = x+e;
    r = r-A*e;
    e = M\r;
    rerr = norm(e)/err0;
    n = n+1;
    if nargout==2
        cv(n+1,:) = [n,x.',rerr];
    end
end
n
if nargout==2
    cv(n+2:end,:) = []; % delete unused space in cv
end
if rerr>TOL && n==nmax
    warning('MYBASIC:noCV','no cv within %i iterations',n)
end

