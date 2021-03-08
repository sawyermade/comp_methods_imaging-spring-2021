alpha = [0 1 2 3];
beta = [-3 : 0.1 : 3];

for i = 1:length(alpha)
    a = alpha(i);
    eigens = zeros([length(beta) 2]);
    
    for j = 1:length(beta)
        b = beta(j);
        M = [a b; b a];
        eigens(j, :) = eig(M);
    end
    
    plot(beta, eigens, '*');
    hold on;
end