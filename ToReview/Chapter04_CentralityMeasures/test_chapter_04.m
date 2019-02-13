%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Masashi Shimbo and Marco Saerens (2018).
% 
% 
% Description: Tests for all code in chapter 04.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format shortE
format compact

% Number of elements for randomply generated tests
n = 10;

rng(72) % set seed of random generator

% % An adjacency matrix for a connected unweighted undirected graph with n nodes
% A = rand(n,n);
% A = A > 0.6325; % sparsify & binarize
% A(1:n+1:end) = 0; % nullify main diagonal (remove self-loops)
% A(2:n+1:end) = 1; % ensure connectedness
% A = double(A | A'); % symmetrize A

% % An adjacency matrix for a connected weighted undirected graph with n nodes
% A = rand(n,n);
% A = (A > 0.6325) .* A; % sparsify
% A(1:n+1:end) = 0; % nullify main diagonal (remove self-loops)
% A(2:n+1:end) = 1; % ensure connectedness
% A = (A + A')/2; % symmetrize A

% Some other test examples

% A = [ 0 1 1 0 0 0
%       1 0 1 0 0 0
%       1 1 0 1 0 0
%       0 0 1 0 1 1
%       0 0 0 1 0 1
%       0 0 0 1 1 0 ];
  
A = [ 0 1 1 1 0 0 0 0 0
      1 0 1 1 0 0 0 0 0
      1 1 0 1 0 0 0 0 0
      1 1 1 0 1 0 0 0 0
      0 0 0 1 0 1 0 0 0
      0 0 0 0 1 0 1 1 1
      0 0 0 0 0 1 0 1 1
      0 0 0 0 0 1 1 0 1
      0 0 0 0 0 1 1 1 0 ];

n = length(A); % number of nodes

% Utilities
e = ones(n, 1); % vector of ones
I = eye(n,n); % identity matrix
H = I - e*e'/n; % centering matrix
eps = 1.0e-100 % precision

sprad = abs(eigs(A, 1)); % spectral radius of A

%% Test 2 Alg_04_02_Brandes(C)

C = 1./(A + eps); % compute cost matrix

addpath('Alg_04_02_Brandes');
bet = Alg_04_02_Brandes(C)
rmpath('Alg_04_02_Brandes');

%% TEST 3 Alg_04_03_Closeness(A, p1, p2)

% hyperparameters
pKatz = 0.5 * (1 / sprad);
pExpm = 1.0;

[cc_k, cc_tc, cc_ks, cc_es] = Alg_04_03_Closeness(A, pKatz, pExpm)

%% TEST 4 Alg_04_04_RandomEccentricity(A)

ec = Alg_04_04_RandomEccentricity(A)

%% TEST 8 Alg_04_08_BagOfPathsBetweenness(A, C, theta)

theta = 2;
C = 1./(A + eps); % compute cost matrix

bopBet = Alg_04_08_BagOfPathsBetweenness(A, C, theta);

%% TEST 9 Alg_04_09_BagOfPathsGroupBetweenness(A, C, hi, hk, theta)

h1 = zeros(n,1); h2 = zeros(n,1);
h1(1) = 1; h1(2) = 1;
h2(7) = 1; h2(8) = 1;

bopGBet = Alg_04_09_BagOfPathsGroupBetweenness(A, C, h1, h2, theta)

%% TEST 10 Alg_04_10_RandomizedShortestPathBetweenness(Pref, C, theta)

theta = 4;
C = 1./(A + eps);  % compute cost matrix
d = A * ones(n,1); % outdegree vector
Dinv = diag(1./d); % outdegree diagonal matrix
Pref = Dinv * A; % transition probabilities matrix of the reference ransom walk

rspBet = Alg_04_10_RandomizedShortestPathBetweenness(Pref, C, theta)

%% TEST 11 Alg_04_11_ApproximateShortestPathBetweenness(A, C, theta)

theta = 20;
C = 1./(A + eps); % compute cost matrix

approxSPBet = Alg_04_11_ApproximateShortestPathBetweenness(A, C, theta)

%% TEST 12 Alg_04_12_RandomizedShortestPathNetFlowBetweenness(Pref, C, theta)

theta = 4;
C = 1./(A + eps);  % compute cost matrix
d = A * ones(n,1); % outdegree vector
Dinv = diag(1./d); % outdegree diagonal matrix
Pref = Dinv * A; % transition probabilities matrix of the reference random walk

rspNetBet  = Alg_04_12_RandomizedShortestPathNetFlowBetweenness(Pref, C, theta)

%% TEST 13 Alg_04_13_CommunicabilityCriticality(A)

cr_communicability = Alg_04_13_CommunicabilityCriticality(A)

%% TEST 14 Alg_04_14_KirchhoffCriticality(A)

cr_Kirchhoff = Alg_04_14_KirchhoffCriticality(A)

%% TEST 15 Alg_04_15_SpanningTreeEdgeCriticality(A)

% % A particular example
% A = [ 0  1  1;
%       1  0  3;
%       1  0  0 ];
  
Cr_st = Alg_04_15_SpanningTreeEdgeCriticality(A)

%% TEST 16 Alg_04_16_BagOfPathsNodeCriticality(A, C, theta)

theta = 2;
C = 1./(A + eps); % compute cost matrix

cr = Alg_04_16_BagOfPathsNodeCriticality(A, C, theta)

%% TEST 17 Alg_04_17_BagOfPathsEdgeCriticality(A, C, theta, Edg)

theta = 3;
C = 1./(A + eps); % compute cost matrix
[row,col,~] = find(A); % compute edge criticality for all edges of G
Edg = [row,col];
Cr = Alg_04_17_BagOfPathsEdgeCriticality(A, C, theta, Edg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
