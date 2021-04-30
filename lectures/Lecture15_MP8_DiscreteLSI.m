%% Lecture 16: Matlab practice 8
clear;
clc;

%% Discrete convolution functions (1D and 2D)

x = randi([0 10], 5,1)
y = randi([-3 4], 3,1)

g = conv(x, y)

% 2D Conv
X = randi([0 10], 5,3)
Y = randi([-3 4], 3,7)

G = conv2(X,Y)


%% 1D LSI system
% h = Delta exp(-Delta^2 (m-n)^2)

Delta = 1.2;
N = 1000
M = 1000
for m=0:M-1
    for n=0:N-1
        h(m+1,n+1) = Delta * exp(-Delta.^2 * (m-n).^2);
    end
end
        
m = 0:1:M-1;
n = 0:1:N-1;

[MM, NN] = meshgrid(m,n);

A = Delta .* exp(-Delta.^2 * (MM-NN).^2);

%  

%% Discretizing 2D LSI system
% Discrete PSF is exp(- \Delta^2 (m1 - n1)^2 - \Delta^2 (m2 - n2)^2)

M = 900;

Delta = .1;
du = 1/(M*Delta);
u = [-M/2:1:(M/2)-1];
[uu, vv] = meshgrid(u, u);

H = exp((vv.^2 + uu.^2) );

% Assuming square pixels
P = sinc(Delta*uu).*sinc(Delta*vv);
 
% Assuming square object discretization
Q = sinc(Delta*uu).*sinc(Delta*vv);


H_sys = (H.*P.*Q);

imagesc(H_sys)





