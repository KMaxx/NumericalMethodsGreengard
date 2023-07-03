function f = exactf(x)
f = -cos(x).*cos(cos(x)) - cos(cos(x)).*sin(x) - sin(cos(x)) - (sin(x)).^2 .* sin(cos(x));