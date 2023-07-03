K1 = 300;
K2 = 1;

m = 20;
t0 = 0;  
u = [3 4 2];
tf = 0.01;  
k = (tf-t0)/(m-1);
t = linspace(0, tf, m);
U(1,1) = 3;
U(2,1) = 4;
U(3,1) = 2;
for n=1:length(t)-1
    rhs1 = u(1)+k*(-K1*u(1)*u(2)+K2*u(3));
    rhs2 = u(2)+k*(-K1*u(1)*u(2)+K2*u(3));
    rhs3 = u(3)+k*(K1*u(1)*u(2)-K2*u(3));
    u = [rhs1; rhs2; rhs3];
    U(1,n+1) = rhs1;
    U(2,n+1) = rhs2;
    U(3,n+1) = rhs3;
end

clf
hold on
plot(t,U(1,1:n+1),'b')
plot(t,U(2,1:n+1),'k')
plot(t,U(3,1:n+1),'r')

legend('u1','u2','u3')
title('u(t) as a function of time')
