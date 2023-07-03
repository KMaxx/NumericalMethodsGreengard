% poisson2.m  -- solve the Poisson problem u_{xx} + u_{yy} = f(x,y)
% on [a,b] x [a,b].  
% 
% The 9-point Laplacian is used at interior grid points.
% This system of equations is then solved using backslash.
% 
% From  http://www.amath.washington.edu/~rjl/fdmbook/chapter3  (2007)


a = 0; 
b = 1; 
m = 40;
h = (b-a)/(m+1);
x = linspace(a,b,m+2);   % grid points x including boundaries
y = linspace(a,b,m+2);   % grid points y including boundaries


[X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
X = X';                     % transpose so that X(i,j),Y(i,j) are
Y = Y';                     % coordinates of (i,j) point

Iint = 2:m+1;              % indices of interior points in x
Jint = 2:m+1;              % indices of interior points in y
Xint = X(Iint,Jint);       % interior points
Yint = Y(Iint,Jint);

f = @(x,y) 1.25*exp(x+y/2)+h^2/12*(1.25*exp(x+y/2)+1.25*exp(x+y/2)/4);         % f_ij function

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
rhs(:,1) = rhs(:,1) - 2*usoln(Iint,1)/(3*h^2)-usoln(Iint-1,1)/(6*h^2)-usoln(Iint+1,1)/(6*h^2);
rhs(:,m) = rhs(:,m) - 2*usoln(Iint,m+2)/(3*h^2)-usoln(Iint-1,m+2)/(6*h^2)-usoln(Iint+1,m+2)/(6*h^2);
rhs(1,:) = rhs(1,:) - 2*usoln(1,Jint)/(3*h^2)-usoln(1,Jint-1)/(6*h^2)-usoln(1,Jint+1)/(6*h^2);
rhs(m,:) = rhs(m,:) - 2*usoln(m+2,Jint)/(3*h^2)-usoln(m+2,Jint-1)/(6*h^2)-usoln(m+2,Jint+1)/(6*h^2);
rhs(1,1) = rhs(1,1) + usoln(1,1)/(6*h^2);
rhs(1,m) = rhs(1,m) + usoln(1,m+2)/(6*h^2);
rhs(m,1) = rhs(m,1) + usoln(m+2,1)/(6*h^2);
rhs(m,m) = rhs(m,m) + usoln(m+2,m+2)/(6*h^2);

% convert the 2d grid function rhs into a column vector for rhs of system:
F = reshape(rhs,m*m,1);

% form matrix A:
I = speye(m);
e = ones(m,1);
T = spdiags([4*e -20*e 4*e],[-1 0 1],m,m);
III = spdiags([e 4*e e],[-1 0 1],m,m);
S = spdiags([e e],[-1 1],m,m);
A = (kron(I,T) + kron(S,III)) / (6*h^2);

% Solve the linear system:
uvec = A\F;  

% reshape vector solution uvec as a grid function and 
% insert this interior solution into usoln for plotting purposes:
% (recall boundary conditions in usoln are already set) 

usoln(Iint,Jint) = reshape(uvec,m,m);

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

axis([a b a b])
daspect([1 1 1])
title('Contour plot of computed solution')
hold off

