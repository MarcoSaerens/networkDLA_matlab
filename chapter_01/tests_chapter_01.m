%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Guillaume Guex (2017).
%
%
% Description: Tests for all code in chapter 01.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test 3

% A simple, connected, directed, weighted graph adjacency matrix with 4 nodes
A = rand(n);
A(A < 0.4) = 0;
 
% the corresponding cost matrix
C = 1 ./ A;

% TEST Alg_01_03_ShortestPathDistance(C)
Alg_01_03_ShortestPathDistance(C)

