function [xk_cgm] = conjugategradientmethod(A,y, num_iter)

[~, Nc] = size(A);

x0 = 0;
AtA = A'*A;
Aty = A'*y;
if nargin==2
num_iter = 100;
end

 
xk_cgm = zeros(Nc,num_iter);
nrk = zeros(num_iter);
err_cgm = zeros(num_iter);
p0 = Aty;
r0 = p0;
 
pk = p0;
rk = r0;


for ii = 1:num_iter
    alphak = norm(rk)^2/(rk'*AtA*pk);
    xk_cgm(:,ii+1) = xk_cgm(:,ii) + alphak.*pk;
    rk1 = rk - alphak*(AtA*pk);
    nrk(ii) =  norm(rk);
    betak1 = norm(rk1)^2/norm(rk)^2;
    pk = rk1 + betak1*pk;
    rk = rk1;
end


end