function sdvscg
format compact
format short e
clc

% 2x2 system
A = [2,-1;-1,2]; b = [0;3];
x = A\b
M = eye(2);
x0 = [1;0];
tol = 1e-5; nmax = 100;
[~,cvsd] = mypsd(A,b,M,x0,tol,nmax);
fprintf('%2i %8.4f %8.4f %9.3e\n',cvsd')

[~,cvcg] = mypcg(A,b,M,x0,tol,nmax);
fprintf('%2i %8.4f %8.4f %9.3e\n',cvcg')

for N = [10,100,1000]
    disp(['N = ' int2str(N)])
    u = ones(N,1);
    A = spdiags([-u,2*u,-u],[-1,0,1],N,N);
    b = randn(N,1);
    tic; x = A\b; toc
    M = speye(size(A));
    x0 = zeros(N,1);
    tol = 1e-10; nmax = 10000000;
    tic; xjac = mybasic(A,b,2*M,x0,tol,nmax); toc
    tic; xsd = mypsd(A,b,2*M,x0,tol,nmax); toc
    tic; xcg = mypcg(A,b,M,x0,tol,nmax); toc
    M = tril(A)*(diag(diag(A))\triu(A)); % SGS prec
    tic; xcg2 = mypcg(A,b,M,x0,tol,nmax); toc
    disp([norm(x-xjac)/norm(x), norm(b-A*xsd)/norm(b)])
    disp([norm(x-xsd)/norm(x), norm(b-A*xsd)/norm(b)])
    disp([norm(x-xcg)/norm(x), norm(b-A*xcg)/norm(b)])
    disp([norm(x-xcg2)/norm(x), norm(b-A*xcg2)/norm(b)])
end

function [x,cv] = mypsd(A,b,M,x,TOL,nmax)
% steepest descent
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
    Ae = A*e;
    alpha = (r'*e)/(e'*Ae);
    x = x+alpha*e;
    r = r-alpha*Ae;
    e = M\r;
    rerr = norm(e,2)/err0;
    n = n+1;
    if nargout==2
        cv(n+1,:) = [n,x.',rerr];
    end
end
if nargout==2
    cv(n+2:end,:) = []; % delete unused space in cv
end
if rerr>TOL && n==nmax
    warning('MYPSD:noCV','no cv within %i iterations',n)
end

function [x,cv] = mypcg(A,b,M,x,TOL,nmax)
% conjugate gradient
x = x(:); % make sure initial guess is column vector
r = b-A*x;
e = M\r;
d = e;
er = e'*r;
err0 = norm(e);  % euclidean magnitude
rerr = 1;        % relative error
n = 0;
if nargout==2
    cv = zeros(nmax+1,2+length(x)); % preallocate space for cv
    cv(1,:) = [n,x.',rerr]; % iter. #, iterate (row vector), rel. error
end
while (rerr>TOL && n<nmax)
    Ad = A*d;
    alpha = er/(d'*Ad);
    x = x+alpha*d;
    r = r-alpha*Ad;
    e = M\r; 
    erold = er; er = e'*r;
    beta = er/erold;
    d = r+beta*d;
    rerr = norm(e,2)/err0;
    n = n+1;
    if nargout==2
        cv(n+1,:) = [n,x.',rerr];
    end
end
if nargout==2
    cv(n+2:end,:) = []; % delete unused space in cv
end
if rerr>TOL && n==nmax
    warning('MYPCG:noCV','no cv within %i iterations',n)
end