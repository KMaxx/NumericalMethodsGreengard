T = 2*pi;
ax = 0;
bx = T;
sigma = 0.7;   % Dirichlet boundary condition at ax
beta = 0.7;     % Dirichlet boundary condtion at bx

m1 = 50;
m2 = m1 + 1;
m = m1 - 1;                 % number of interior grid points
h = (bx-ax)/m1;  % average grid spacing, for convergence tests

% set grid points:  
gridchoice = 'uniform';          % see xgrid.m for other choices
% gridchoice = 'random'; 
t = xgrid(ax,bx,m,gridchoice);   

% set up matrix A (using sparse matrix storage):
A = spalloc(m2,m2,3*m2);   % initialize to zero matrix

% first row for Dirichlet BC at ax:
A(1,1:3) = fdcoeffF(0, t(1), t(1:3)); 

% interior rows:
for i=2:m1
    A(i,i-1:i+1) = fdcoeffF(2, t(i), t((i-1):(i+1)));
end

% last row for Dirichlet BC at bx:
A(m2,m:m2) = fdcoeffF(0,t(m2),t(m:m2)); 
G = @(x) A*x + sin(x);  % true soln

% set up matrix JM (using sparse matrix storage):
JM = spalloc(m2,m2,3*m2);   % initialize to zero matrix

JM(1,1:3) = fdcoeffF(2, t(1), t(1:3)); 

% interior rows:
for i=2:m1
    JM(i,i-1:i+1) = fdcoeffF(2, t(i), t((i-1):(i+1)));
end

JM(m2,m:m2) = fdcoeffF(2,t(m2),t(m:m2)); 
J = @(x) JM + diag(cos(x));  % true soln

theta = 0.7*cos(t);
% theta = 0.7 + sin(t/2);
% theta = 0.7*cos(t) + 0.5*sin(t);
% theta = 0.7*ones(size(t));
hold on
plot(t,theta)  % plot true solution
hold off
for k=1:10
    % solve linear system:
    J1 = J(theta);
    J1 = J1(2:m2-1,2:m2-1);
    G1 = G(theta);
    G1 = G1(2:m2-1);
    G1 = -G1;
    dtheta = J1\G1;
    theta(2:m2-1) = theta(2:m2-1) + dtheta;
    hold on
    plot(t,theta)  % plot true solution
    hold off
end
