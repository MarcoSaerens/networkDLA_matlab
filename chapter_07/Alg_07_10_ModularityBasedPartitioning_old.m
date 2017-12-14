function v1 = Alg_07_10_ModularityBasedPartitioning_old(A)
%   Modularity-based two-way partitioning of a graph
%       Chapter 7 : Clustering nodes
%           Section 9 : Modularity criterion and its spectral optimization
%               Subsection 3 : Two-way partitioning based on modularity
%
% Arguments:
%
% - A: the n x n symmetric adjacency matrix, associated with the
%      weighted undirected graph G containing n nodes and no dangling node
%
% Returns:
%
% - v1: the n x 1 leading eigenvector (correponding to the largest 
%       eigenvalue) of the modularity matrix Q
%
% (c) Simon Gilissen & Edouard Lardinois

% Compute the degree vector
d = sum(A,2);

% Compute the volume of G
vol_G = sum(d);

% Compute the modularity matrix
Q = A - d*d'/vol_G;

% Compute the eigenvectors of Q : Q v = lambda v
[EV,D] = eig(Q);
lambda = diag(D);
% Sort the eigenvectors and eigenvalues in decreasing order of eigenvalue
% and keep the dominant one.
[~, lambda_idx] = sort(lambda,'descend');
EV = EV(:,lambda_idx);
v1 = EV(:,1);

