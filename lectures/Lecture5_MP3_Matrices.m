%% Lecture 4: MP3 - Manipulating matrices
%
%% 1. Creating matrices and indexing
A = [ 1 -1 9; -7 0 4]

entry22 = A(2,2)

%% 2. Matrix addition and subtraction

B = [7 7 7; 7 7 7];
B1 = 7*ones(2,3);

AplusB = A + B;

BminusA = B - A;

%% 3. Manipulying matrices by:
%   a) Scalars
%   b) Vectors
%   c) Matrices (Matlab will throw an error, for dimension issues)

MatrixB = [1 -4 4; 9 -3 -2; 7 7 7; 7 7 7] 

pitimesMatrixB = pi * MatrixB;

Vectora = [2 ; 3; -1 ];

matVecMult = MatrixB * Vectora

Vectorb = [1 ; 1];

%matVecMult = MatrixB * Vectorb % Not possible

MatrixB
A'

matMatMult = MatrixB * (A')

elementwiseprod = A .* B


A.^2

%% 4. Computing inner products and norms of vectors
clc;
norm(Vectora)

Vectorb = [ 1 1 1]';
sqrt(Vectora' * Vectora)

Vectora' * Vectorb

%% 5. Checking linear independence
%    a) Rank and determinant of matrices
%    b) Can also look at reduced row echelon form


