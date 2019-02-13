%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Guillaume Guex (2017).
%
%
% Description: Tests for all code in Chapter 02.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of nodes for tests
n = 10;

%% Test 1

% A simple, connected, undirected, unweighted graph adjacency matrix with
% n nodes
A = rand(n) > 0.4;
A(1:n+1:end) = 0;
A = (A + A') > 0;

% TEST Alg_02_01_LocalSimilarityMeasures(A, i, k)
sim = Alg_02_01_LocalSimilarityMeasures(A, 1, 4);
sim.direct
sim.common
sim.pref
sim.cos
sim.Jaccard
sim.Dice
sim.prohub
sim.dehub

%% Test 2

% A simple, connected, undirected, weighted graph adjacency matrix with
% n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;
A = (A + A')/2;

% TEST Alg_02_02_KatzSimilarityAndLeichtsExtension(A, alpha)
K = Alg_02_02_KatzSimilarityAndLeichtsExtension(A, rand(1)/10);
K.Katz
K.Leicht

%% Test 3

% A simple, connected, undirected, weighted graph adjacency matrix with
% n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;
A = (A + A')/2;

% TEST Alg_02_03_CommuteTimeAndEuclideanCommuteTimeDistances(A)
D_com = Alg_02_03_CommuteTimeAndEuclideanCommuteTimeDistances(A);
D_com.CT
D_com.ECT
D_com.CCT
D_com.CECT


%% Test 4

% A simple, connected, directed, weighted graph adjacency matrix with
% n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;

% TEST Alg_02_04_SimRankSimilarity(A, alpha)
Alg_02_04_SimRankSimilarity(A, rand(1))

%% Test 5

% Two simple, connected, directed, weighted graph adjacency matrices
% with n and n+2 nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;
B = rand(n+2);
B(1:n+3:end) = 0;
B(B < 0.4) = 0;

% TEST Alg_02_05_BlondelSimilarity(A, B)
Alg_02_05_BlondelSimilarity(A, B)

%% Test 6

% A simple, connected, undirected, weighted graph adjacency matrix with
% n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;
A = (A + A')/2;

% TEST Alg_02_06_ExponentialDiffusionKernel(A, alpha)
K_exp = Alg_02_06_ExponentialDiffusionKernel(A, rand(1) / 10 );
K_exp.ED
K_exp.LED

%% Test 7

% A simple, connected, undirected, weighted graph adjacency matrix with
% n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;
A = (A + A')/2;

% TEST Alg_02_07_ModifiedRegularizedLaplacianKernel(A, alpha, gamma)
Alg_02_07_ModifiedRegularizedLaplacianKernel(A, rand(1), rand(1))

%% Test 8

% A simple, connected, undirected, weighted graph adjacency matrix with
% n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;
A = (A + A')/2;

% TEST Alg_02_08_CommuteTimeKernel(A)
K_com = Alg_02_08_CommuteTimeKernel(A);
K_com.CT
K_com.CCT

%% Test 9

% A simple, connected, undirected, weighted graph adjacency matrix with
% n nodes
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;
A = (A + A')/2;

% TEST Alg_02_09_RegularizedCommuteTimeKernel(A, alpha)
K_rctwr = Alg_02_09_RegularizedCommuteTimeKernel(A, rand(1));
K_rctwr.RCT
K_rctwr.RWR

%% Test 10

% A simple, connected, directed, weighted graph adjacency matrix with
% n nodes and a weight vector
A = rand(n);
A(1:n+1:end) = 0;
A(A < 0.4) = 0;
w = rand(n,1);
w = w / sum(w);

% TEST Alg_02_10_MarkovDiffusionSquareDistance(A, w, t)
st = Alg_02_10_MarkovDiffusionSquareDistance(A, w, 4);
st.D_MD
st.D_MDA
st.K_MD
st.K_MDA
