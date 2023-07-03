maxiter = 100;        % number of iterations to take

m = 199;
ax = 0;
bx = 1;
alpha = 0;
beta = 0;
f = @(x) ones(size(x));   % f(x) = 1

h = (bx-ax) / (m+1);
x = linspace(ax, bx, m+2)';

% determine exact solution to linear system by setting up
% system Au = f and solving with backslash:
e = ones(m+2,1);

S = 1/h^2 * spdiags([e -2*e e], [-1 0 1], m+2, m+2);
S(1,1:2) = [1 0];
S(m+2,m+1:m+2) = [0 1];
F = 1/(2*h) * spdiags([-e e], [-1 1], m+2, m+2);
F(1,1:2) = [0 0];
F(m+2,m+1:m+2) = [0 0];
I = speye(m+2);
I(1,1:2) = [0 0];
I(m+2,m+1:m+2) = [0 0];
A = -S + F + I;
%[L,U] = ilu(A);
rhs = f(x);
rhs(1) = alpha;
rhs(m+2) = beta;
ustar = A\rhs;
restart1 = 20;
restart2 = [];
tol = [];
maxit1 = 5;
maxit2 = 100;
[u,flag1,relres,iter,resvec1]=gmres(A,rhs,restart1,tol,maxit1);
[u,flag2,relres,iter,resvec2]=gmres(A,rhs,restart2,tol,maxit2);
figure(1)
clf
plot(resvec1,'r')
hold on
plot(resvec2,'b')
title('Relative Residuals','FontSize',15)
xlabel('iteration')
ylabel('relative residual')
legend('with restarting','without restarting')