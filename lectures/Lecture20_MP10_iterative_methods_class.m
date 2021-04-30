%% Lecture20_MP10_iterative_methods
%% Fixed point iterations: Example using a highly non-linear equation
clear;
clc;

x = linspace(0,3, 100);

fx = x;
gx = @(x)( exp(-x.^2)./sin(x+1) );
num_iter = 20;

figure(1); clf
plot(x,fx,'linewidth', 3); hold on;
plot(x,gx(x),'linewidth', 3);

x1 = 2;
x_est(1) = x1;
for k=2:num_iter
    x_est(k) = gx(x_est(k-1));
    
    
    line([x_est(k-1) x_est(k-1)],[x_est(k-1) x_est(k)],'color','k');
    plot(x_est(k-1),x_est(k),'.k','markersize',20,'linewidth',20); 
    pause(1); drawnow;
    line([x_est(k-1) x_est(k)],[x_est(k) x_est(k)],'color','k');
    pause(1); drawnow;
end

% x_est(90:100)




%% Steepest descent

Nr = 40;
Nc = 25;
num_iter = 50;

A = randn(Nr, Nc);
y = randi([-20 20], Nr, 1);

AtA = A'*A;
Aty = A'*y;

x_est_k = zeros(Nc, num_iter);
xtrue = A\y;

for kk=1:num_iter
    
    rk = AtA*x_est_k(:,kk) - Aty;
    nrk(kk) = norm(rk).^2;
    nArk(kk) = norm(A*rk).^2;
    lambda_k = nrk(kk)/nArk(kk);
    x_est_k(:, kk+1) = x_est_k(:,kk) - lambda_k.*rk;
    error_sd(kk) = norm(xtrue-x_est_k(:,kk+1));
end

figure(2); clf;
plot(x_est_k(:,end));
hold on;
plot(A\y,'o');




 
 
%% Conjugate gradient method

[xk_cgm] = conjugategradientmethod(A, y, num_iter);
error_cgm= sqrt(sum((xtrue*ones(1, num_iter + 1)-xk_cgm).^2,1));

 
figure(3);
plot(xk_cgm(:,end), 'r', 'linewidth', 2.5)
 
figure(4); clf;
plot((error_sd), 'linewidth', 2.5); hold on;
plot((error_cgm), 'linewidth', 2.5);



