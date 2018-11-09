%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Guillaume Guex (2017).
%
%
% Description: Tests for all code in chapter 06.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of nodes for tests
n = 5;

%% Test 1

% A simple, connected, undirected, weighted graph adjacency matrix with 
% n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.3) = 0;
A = (A + A')/2;

Y = zeros(n, 2);
Y(2, 1) = 1;
Y(4, 2) = 1;

% TEST Alg_06_01_SimpleRegularizationClassification(A, Y, lambda)
Alg_06_01_SimpleRegularizationClassification(A, Y, 0.5)

%% Test 2

Y = zeros(n, 3);
Y(1, 1) = 1;
Y(2, 2) = 1;
Y(3, 3) = 1;

% TEST Alg_06_02_HarmonicFunctionClassification(A, Y)
Alg_06_02_HarmonicFunctionClassification(A, Y)

%% Test 3

Y = zeros(n, 2);
Y(1, 1) = 1;
Y(2, 1) = 1;
Y(4, 2) = 1;

% TEST Alg_06_03_RWWRClassification(A, Y, alpha)
Alg_06_03_RWWRClassification(A, Y, 0.9)

%% Test 4

% TEST Alg_06_04_DWalkClassification(A, Y, alpha)
Alg_06_04_DWalkClassification(A, Y, 0.6)

%% Test 5

Y = zeros(n, 2);
Y(1, 1) = 1;
Y(2, 1) = 1;
Y(4, 2) = 1;
Y(5, 2) = 1;
C = 1 ./ A;

% TEST Alg_06_05_BoPBetweennessClassification(A, C, Y, theta)
Alg_06_05_BoPBetweennessClassification(A, C, Y, 0.1)

%% Test 6
X = rand(n, 3);
Y = zeros(n, 2);
Y(2, 1) = 1;
Y(4, 2) = 1;

% TEST Alg_06_06_RidgeWithLaplacianRegClassification(A, X, Y, lambda, mu)
Alg_06_06_RidgeWithLaplacianRegClassification(A, X, Y, 0.4, 0.2)

%% Test 7
K = (eye(n) - ones(n,n)/n)*X*X'*(eye(n) - ones(n,n)/n);

% TEST Alg_06_07_KernelRidgeWithLaplacianRegClassification(A, K, Y, lambda, mu)
Alg_06_07_KernelRidgeWithLaplacianRegClassification(A, K, Y, 0.4, 0.2)
