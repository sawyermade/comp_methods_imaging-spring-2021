%% Lecture 6 (completed 7): MP3 - Manipulating matrices
% Matrices
%% Matrix transpose (adjoint)

A = randi([-20 20], 4, 4)
Atranspose = A'

%% Matrix inverse (square matrix)
rank(A)
det(A)
Ainv = inv(A);

Ainv*A % left inv

A*Ainv % right inverse


%% Left inverses (rectangular matrix, tall matrix)

B = randi([-20 20], 10, 4)

Bleftinv = pinv(B) % pseudo-inverse

Bleftinv * B


B*Bleftinv

%% Right inverses (fat matrix)

C = randi([-20 20], 10, 20);
Crightinv = pinv(C)

C*Crightinv 




%% Eigenvalue decomposition
% Eigenvalues & Eigenvectors
%
% 

[V, S] = eig(A)


%% Plotting 2D functions

x= -10:1:10
y= -10:1:10



[XX, YY] = meshgrid(x,y);

fxy = exp(-(XX.^2 + YY.^2)/20);

surf(XX,YY, fxy);
% shading interp





















