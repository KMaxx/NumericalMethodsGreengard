function [f,fx,err] = fftpdesolver2(a,b,N,alpha, betha)
x = sin((a-(b-a)/2:(b-a)/N:b-(b-a)/2-(b-a)/N)).^2;
y = fft(x);
index = (2*pi/(b-a))*[-N/2:N/2-1];
index = fftshift(index);
yx = i*index.*y;
yxx = -index.*index.*y;
f = real(ifft(yxx + alpha*yx - betha*y));
fx = exactf2((a-(b-a)/2:(b-a)/N:b-(b-a)/2-(b-a)/N));
err = immse(f,fx);
plot((a-(b-a)/2:(b-a)/N:b-(b-a)/2-(b-a)/N), f);
%plot((a-(b-a)/2:(b-a)/N:b-(b-a)/2-(b-a)/N), fx);
