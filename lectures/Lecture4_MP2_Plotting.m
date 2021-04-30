%% Lecture 4: MP2 - Plotting 1D functions
%
%% 1. Plot x(t) = t
t = [-1:0.1:1]';

mylinearfunction = t;

h1 = plot(t, mylinearfunction);

hold on;
%% 2. Plot x(t) = cos(pi*t)

t = [-1:0.05:1];
oursecondfunction = cos( pi*t );
plot(t, oursecondfunction);
grid on;
