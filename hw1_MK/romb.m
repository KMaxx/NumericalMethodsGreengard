function [I, T, R, MI, ES, E] = romb(N,a,b,f)
T = zeros(N,N);
R = zeros(N,N);
MI = zeros(N,N);
ES = zeros(N,N);
E = ones(N,N)*integral(f,a,b);
for m = 1:N
    M = 2^(m-1);
    h = (b-a)/M;
    x = f(a:h:b);
    x(1) = x(1)/2;
    x(end) = x(end)/2;
    T(m,1) = h*(sum(x));
    MI(m,1) = M;
    ES(m,1) = 1;
    for n = 2:m
        T(m,n) = T(m,n-1) + (T(m,n-1)-T(m-1,n-1))/(4^(n-1)-1);
        ES(m,n) = ES(m,n-1)+ES(m-1,n-1);
        MI(m,n) = M;
    end
end 
E = E-T;
for m = 1:N-1
    for n = 1:m-1
        R(m,n) = (T(m,n)-T(m-1,n))/(T(m+1,n)-T(m,n));
    end
end

I = T(N,N);