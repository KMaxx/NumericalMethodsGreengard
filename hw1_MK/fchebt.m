function [u,err] = fchebt(N)
t = -cos(((0:N-1)+1/2)/N*pi);
ft = cos(t).^2;
alpha = zeros([1 N]);
for n = 1:N
    alpha(n) = 2/N*dot(ft,cos((n-1)*acos(t)));
end
betha = zeros([1 N]);
c = 0.0;
sign = -1;
for n = 2:N-1
    betha(n) = (alpha(n-1)-alpha(n+1))/2/(n-1);
    c = c - sign*betha(n);
    sign = (-1)*sign;
end
betha(N) = alpha(N-1)/2/(N-1);
c = c - sign*betha(N);
betha(1) = c;
u = zeros([1 N]);

for n = 1:N
    u(n) = dot(betha,cos((0:N-1)*acos(t(n))))+1;
end
err = immse(u,(t+sin(2*t)/2)/2+1-(-1+sin(-2)/2)/2);
plot(t, u)