function H = H_system(D, d, zed)
    M = 1000;
    N = 1000;
    lambda = 0.5 * 10^-6;
    z = zed * 10^-3;
    
    Delta = D * 10^-6;
    delta = d * 10^-6;
    
    dm = 1/(M * delta)
    dn = 1/(N * Delta)
    
    m = [-M/2:1:(M/2)-1] .* dm;
    n = [-N/2:1:(N/2)-1] .* dn;
    
    [mm, nn] = meshgrid(m, n);
    
    H = exp((-1 * sqrt(-1) * pi * lambda * z) .* (mm.^2 + nn.^2));
    
    P = sinc(Delta * mm) .* sinc(delta * mm)
    
    H_sys = H .* P .* Q;
    
    figure;
    imagesc(real(H_sys));
    
    figure;
    surf(abs(real(H_sys)), 'edgecolor', 'none');
end