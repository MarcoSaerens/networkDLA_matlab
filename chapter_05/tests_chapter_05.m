%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: François Fouss (2017).
%
%
% Description: Tests for all code in chapter 05.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% COMPUTE Alg_01_03_ShortestPathDistance(C)
D = Alg_01_03_ShortestPathDistance(C);
% TEST Alg_05_01_ProximityPrestigeScores(D)
Alg_05_01_ProximityPrestigeScores(D)