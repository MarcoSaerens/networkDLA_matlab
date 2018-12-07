%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Masashi Shimbo (2018).
% 
% 
% Description: Tests for all code in chapter 04.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of elements for tests
n = 10;

% utilities
e = ones(n, 1); % vector of ones
I = eye(n); % identity matrix
H = I - e*e'/n; % centering matrix
eps = 10^(-100); % precision

rng(72) % set seed of random generator

% An adjacency matrix for a connected undirected graph with n nodes
A = rand(n);
A = A > 0.6325; % sparsify & binarize
A(1:n+1:end) = 0; % remove diagonal (self-loops)
A = real(A | A'); % symmetrize

% A
sprad = eigs(A, 1); % spectral radius of A

% Test 1 Alg_04_02_Brandes(C)

C = 1./(A + eps); % compute cost matrix

addpath('Alg_04_02_Brandes');
bet = Alg_04_02_Brandes(C)
rmpath('Alg_04_02_Brandes');

% TEST 2 Alg_04_03_Closeness(A, p1, p2)

% hyperparameters
pKatz = 0.5 * (1 / sprad);
pExpm = 1.0;

[cc_k, cc_tc, cc_ks, cc_es] = Alg_04_03_Closeness(A, pKatz, pExpm)

% TEST 3 Alg_04_04_RandomEccentricity(A)
ec = Alg_04_04_RandomEccentricity(A)

% TEST 8 Alg_04_08_BagOfPathsBetweenness(A, C, theta)

theta = 2;
C = 1./(A + eps); % compute cost matrix

bopBet = Alg_04_08_BagOfPathsBetweenness(A, C, theta)

% TEST 9 Alg_04_09_BagOfPathsGroupBetweenness(A, C, hi, hk, theta)

h1 = zeros(n,1); h2 = zeros(n,1);
h1(1) = 1; h1(2) = 1;
h2(7) = 1; h2(8) = 1;

bopGBet = Alg_04_09_BagOfPathsGroupBetweenness(A, C, h1, h2, theta)

