%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: François Fouss (2017).
%
%
% Description: Tests for all code in chapter 05.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of nodes for tests
n = 10;

% %% Test 1
% 
% % A simple, connected, directed, weighted adjacency matrix with n nodes
% A = rand(n);
% A(1:n+1:end) = 0;
% A(A < 0.4) = 0;

%% Test 1

% A simple, connected, directed, weighted adjacency matrix with n nodes
n = 10;
A = rand(n);
A(A < 0.4) = 0;

%A =[0 0 0 0 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0 0 0; 1 1 0 0 0 0 0 0 0 0 0; 0 1 0 1 0 1 0 0 0 0 0; 0 1 0 0 1 0 0 0 0 0 0; 0 1 0 0 1 0 0 0 0 0 0; 0 1 0 0 1 0 0 0 0 0 0; 0 1 0 0 1 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0];
%A = [0 1 0 0 0 1; 0 0 1 1 0 0; 0 0 0 1 1 1; 1 0 0 0 0 0; 0 0 0 0 0 0; 1 0 0 0 0 0];

%A = [0 1 1 1; 0 0 1 1; 0 1 0 0; 0 0 0 0]

% the corresponding cost matrix
C = 1 ./ A;

% Compute Alg_01_03_ShortestPathDistance(C)
D = Alg_01_03_ShortestPathDistance(C);

% TEST Alg_05_01_ProximityPrestigeScores(D)
Alg_05_01_ProximityPrestigeScores(D)

%% Test 2

% A simple, connected, directed, weighted adjacency matrix with n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;

% TEST Alg_05_02_Bonacich(A)
Alg_05_02_Bonacich(A)

%% Test 3

% A simple, connected, directed, weighted adjacency matrix with n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;

u = ones(n,1)/n;

% TEST Alg_05_03_KatzHubbell(A,u,alpha)
Alg_05_03_KatzHubbell(A,u,0.1)

%% Test 4

% A simple, connected, directed, weighted adjacency matrix with n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;

% TEST Alg_05_04_CitationInfluence(A)
Alg_05_04_CitationInfluence(A)

%% Test 5

% A simple, connected, directed, weighted adjacency matrix with 4 nodes
A = [0 4 6 0; 1 0 0 2; 2 0 0 0; 0 3 0 0];
eps = [1 2; 2 1; 1 3; 3 1; 2 4; 4 2];

% TEST Alg_05_05_LeastSquaresRating(A,eps)
Alg_05_05_LeastSquaresRating(A,eps)

%% Test 6

% A simple, connected, directed, weighted adjacency matrix with n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;

% Dangling nodes are made absorbing
d  = sum(A,2); % Outdegree vector
if min(d) == 0
    A(d,d) = 1;
end

u = ones(n,1)/n;

% TEST Alg_05_06_PageRank(A, u, alpha)
Alg_05_06_PageRank(A, u, 0.85)

%% Test 7

% A simple, connected, directed, weighted adjacency matrix with n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;

% TEST Alg_05_07_HITS(A)
Alg_05_07_HITS(A)
