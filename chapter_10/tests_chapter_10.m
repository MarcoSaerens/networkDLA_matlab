%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Guillaume Guex (2017).
% 
% 
% Description: Tests for all code in chapter 10.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of elements for tests
n = 5;

% utilities
e = ones(n, 1);
I = eye(n);
H = I - e*e'/n;

%% Test 1

% Coordinate matrix
X_init = rand(n, 2);
% Covariance matrix 
K_init = H*X_init*X_init'*H;

% TEST Alg_10_01_KernelPrincipalComponentAnalysis(K, p)
Alg_10_01_KernelPrincipalComponentAnalysis(K_init, 2)

%% Test 2

% Euclidean distance matrix 
D_init = sqrt( diag(K_init)*e' + e*diag(K_init)' - 2*K_init );

% TEST Alg_10_02_ClassicalMultidimensionalScaling(D, p)
Alg_10_02_ClassicalMultidimensionalScaling(D_init, 2)

%% Test 3

% A simple, connected, undirected, weighted graph adjacency matrix with 
% n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;
A = (A + A')/2;

% TEST Alg_10_03_MarkovDiffusionMap(A, t, p)
Alg_10_03_MarkovDiffusionMap(A, 3, 2)

%% Test 4

% TEST Alg_10_04_LaplacianEigenmapEmbedding(A, p)
Alg_10_04_LaplacianEigenmapEmbedding(A, 2)

%% Test 5

% transform A to binary matrix
A_bin = A > 0;

% compute the MDS embedding with SP-distance
D_sp = 1 ./ A_bin;
D_sp(1:n+1:end) = 0;
for k = 1:n
     D_sp = min( D_sp, (D_sp(:, k)*(e') + e*D_sp(k, :)) );
end
K_sp = -1/2 .* H*D_sp*H;
[U, Lambda] = eig(K_sp);
X_0 = real(U*sqrt(Lambda));

% TEST Alg_10_05_LatentSocialMap(A, p)
Alg_10_05_LatentSocialMap(A_bin, n, X_0)

%% Test 6

% TEST Alg_10_06_SpringNetworkLayout(Delta, X_0, l_0, k)
Alg_10_06_SpringNetworkLayout(D_sp, X_0, 1, 1)

%% Test 7

% weight vector and matrix
W = rand(n) .* A_bin;
W = W + W';
w = rand(n, 1);

% TEST Alg_10_07_ForcedirectedLayoutGraph(W, w, X_0, a, r)
Alg_10_07_ForcedirectedLayoutGraph(W, w, X_0, 2, -1)
