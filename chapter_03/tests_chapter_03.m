%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Guillaume Guex (2017).
%
%
% Description: Tests for all code in chapter 03.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of nodes for tests
n = 5;

%% Test 1

% A simple, connected, undirected, weighted graph adjacency matrix with 
% n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;
A = (A + A')/2;

% TEST Alg_03_01_LogarithmicForestDistance(A, alpha)
Alg_03_01_LogarithmicForestDistance(A, 0.1)

%% Test 2

% A simple, connected, directed, weighted graph adjacency matrix with 
% n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.2) = 0;

% The corresponding cost matrix
C = 1 ./ A;

% TEST Alg_03_02_RegularBagOfPathsProbability(A, C, theta)
Alg_03_02_RegularBagOfPathsProbability(A, C, 0.1)

%% Test 3

% TEST Alg_03_03_BagOfHittingPathsProbability(A, C, theta)
Pi_h = Alg_03_03_BagOfHittingPathsProbability(A, C, 0.1);
Pi_h.std
Pi_h.bar

%% Test 4

% TEST Alg_03_04_BagOfHittingPathsSurprisalDistance(Pi_h)
Alg_03_04_BagOfHittingPathsSurprisalDistance(Pi_h.std)

%% Test 5

% TEST Alg_03_05_BagOfHittingPathsFreeEnergyDistance(A, theta, C)
Alg_03_05_BagOfHittingPathsFreeEnergyDistance(A, C, 0.1)

%% Test 6

% TEST Alg_03_06_RandomizedShortestPathDissimilarity(A, theta, C)
Alg_03_06_RandomizedShortestPathDissimilarity(A, C, 0.1)

%% Test 7

% TEST Alg_03_07_BagOfPathsAbsorptionProbability(A, C, theta, abs)
Alg_03_07_BagOfPathsAbsorptionProbability(A, C, 0.1, [0, 1, 0, 1, 0]')

%% Test 8
P_ref = diag(A * ones(n, 1))^(-1) * A;

% TEST Alg_03_08_BagOfPathsCovariance(C, P_ref, theta)
Alg_03_08_BagOfPathsCovariance(C, P_ref, 0.1)