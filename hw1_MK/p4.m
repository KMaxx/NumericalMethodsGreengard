N = 4:2:16;
err = zeros(length(N));
for n = 1:length(N)
    [u,er] = fchebt(N(n));
    err(n) = log(vpa(er));
end
plot(N,vpa(err));