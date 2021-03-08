alpha = [0 : 0.1 : 5];
beta = [0 : 0.1 : 5];

[a, b] = meshgrid(alpha, beta);
abdet = a*a - b*b;
surf(a, b, abdet);