c1 = 2;
ax = 0;
bx = pi;

for c2=-10:2:10
    utrue = @(x) c2*sin(x) + c1*cos(x);  % true soln

    % true solution on fine grid for plotting:
    xfine = linspace(ax,bx,101);
    ufine = utrue(xfine);
    hold on
    plot(xfine,ufine)  % plot true solution
    hold off
end