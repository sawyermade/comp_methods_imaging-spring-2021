%% Manipulating vectors
%{
1. Creating vectors

2. Manipulating vectors
    a) Scalar multiplication
    b) Vector addition
    c) Indexing

3. Computing inner products and norms of vectors

4. Checking linear independence
    a) Rank and determinant of matrix representation of the vectors
    b) Can also look at reduced row echelon form
%}

%%
clc;

% 1)
v = [0; 3; 2];

%2(a) Scalar mult

u = 5*v;

% 2(b)

w = [1 2 3]';

v;

vplusw = v + w;

% 2(c)

w(1) ;

% 3:

innerprod_vw = v' * w;

innerprod2_vw = dot(v, w);


% 4:
v1 = [1 -2 0]'

v2 = [4 0 8]'

v3 = [3 -1 5]'


A = [v1 v2 v3]

determinantA = det(A)

matB = [1 0 0; 0 1 0; 1 1 1]

det(matB)

