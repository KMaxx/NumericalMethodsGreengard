% poisson1_2.m  -- solve the Poisson problem u_{xx} + u_{yy} = f(x,y)
% on [a_x,b_x] x [a_y,b_y] with dx=dy=h.  
% 
% The 5-point Laplacian is used at interior grid points.
% This system of equations is then solved using backslash.
% 
% From  http://www.amath.washington.edu/~rjl/fdmbook/chapter3  (2007)


a = 0; 
b = 1;
ay = 1;
by = 4;
m = 20;
h = (b-a)/(m+1);
n = round((by-ay)/h-1);  % choose n such that dx=dy
x = linspace(a,b,m+2);   % grid points x including boundaries
y = linspace(ay,by,n+2);   % grid points y including boundaries


[X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
X = X';                     % transpose so that X(i,j),Y(i,j) are
Y = Y';                     % coordinates of (i,j) point

Iint = 2:m+1;              % indices of interior points in x
Jint = 2:n+1;              % indices of interior points in y
Xint = X(Iint,Jint);       % interior points
Yint = Y(Iint,Jint);

f = @(x,y) 1.25*exp(x+y/2);         % f(x,y) function

rhs = f(Xint,Yint);        % evaluate f at interior points for right hand side
                           % rhs is modified below for boundary conditions.

utrue = exp(X+Y/2);        % true solution for test problem

% set boundary conditions around edges of usoln array:

usoln = utrue;              % use true solution for this test problem
                            % This sets full array, but only boundary values
                            % are used below.  For a problem where utrue
                            % is not known, would have to set each edge of
                            % usoln to the desired Dirichlet boundary values.


% adjust the rhs to include boundary terms:
rhs(:,1) = rhs(:,1) - usoln(Iint,1)/h^2;
rhs(:,n) = rhs(:,n) - usoln(Iint,n+2)/h^2;
rhs(1,:) = rhs(1,:) - usoln(1,Jint)/h^2;
rhs(m,:) = rhs(m,:) - usoln(m+2,Jint)/h^2;


% convert the 2d grid function rhs into a column vector for rhs of system:
F = reshape(rhs,m*n,1);

% form matrix A:
I = speye(m);
In = speye(n);
e = ones(m,1);
en = ones(n,1);
T = spdiags([e -4*e e],[-1 0 1],m,m);
S = spdiags([en en],[-1 1],n,n);
A = (kron(In,T) + kron(S,I)) / h^2;


% Solve the linear system:
uvec = A\F;  

% reshape vector solution uvec as a grid function and 
% insert this interior solution into usoln for plotting purposes:
% (recall boundary conditions in usoln are already set) 

usoln(Iint,Jint) = reshape(uvec,m,n);

% assuming true solution is known and stored in utrue:
err = max(max(abs(usoln-utrue)));   
fprintf('Error relative to true solution of PDE = %10.3e \n',err)

% plot results:

clf
hold on

% plot grid:
% plot(X,Y,'g');  plot(X',Y','g')

% plot solution:
contour(X,Y,usoln,30,'k')

axis([a b ay by])
daspect([1 (by-ay)/(b-a) 1])
title('Contour plot of computed solution')
hold off

