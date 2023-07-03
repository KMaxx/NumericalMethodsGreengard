N = 2:2:20;
err = zeros(length(N));
for n = 1:length(N)
    [f,fx,er] = fftpdesolver2(0,2*pi, N(n), 1, 1);
    err(n) = log(vpa(er));
end
plot(N,vpa(err));