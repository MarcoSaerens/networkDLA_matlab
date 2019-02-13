%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Guillaume Guex (2017).
%
%
% Description: Tests for all code in Chapter 01.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test 2-3

% A connected, directed, weighted adjacency matrix with n nodes
n = 10;
A = rand(n);
A(A < 0.4) = 0;
 
% The corresponding cost matrix
C = 1 ./ A;

% TEST Alg_01_02_ShortestPathDistanceElementwise(C)
Alg_01_02_ShortestPathDistanceElementwise(C)

% TEST Alg_01_03_ShortestPathDistance(C)
Alg_01_03_ShortestPathDistance(C)

