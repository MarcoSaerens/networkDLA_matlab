function [ bet ] = Alg_04_08_BagOfPathsBetweenness(A, C, theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Bertrand Lebichot, revised by Marco Saerens (2018).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016),
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: The bag of paths approach for computing a betweenness
% centrality measure on the nodes of a weighted directed or undirected,
% strongly connected, graph G without self-loops.
%
% INPUT:
% ------
% - A: the n x n adjacency matrix of a weighted directed or undirected,
%   strongly connected, graph G containing n nodes.
% - C: the n x n cost matrix C associated to G (if not specified,
%   the costs are the inverse of the affinities, but other 
%   choices are possible).
% - theta: the (scalar) inverse temperature parameter.
%
% OUTPUT:
% ------- 
% - bet: the nx 1 bag of hitting paths betweenness vector  
%   containing the probabilities that a randomly chosen path 
%   visits an intermediate node.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check if A is a square matrix 
[n,m] = size(A);
if n ~= m
    display('Error : matrix A must be square');
    return;
end

I  = eye(n); % identity matrix (same size as A)
%% Algorithm
d = A * ones(n,1); % the outdegree vector
Dinv = diag(1./d); % the inverse outdegree square matrix
Pref = Dinv * A;  % the reference transition probabilities matrix

W = Pref .* exp(-theta*C); % the W matrix

Z  = (I - W)\I; % the fundamental matrix Z
Z0 = Z - diag(diag(Z)); % set diagonal to zero

Dzinv = diag(diag(Z).^-1); % the inverse of diagonal matrix Dz

N = Z0 * Dzinv * Z0; % matrix of normalization factors
Ndiv = 1./(N + eps); % matrix T contains elements 1/nij

% Computation of the BoP betweenness
bet = Dzinv * diag( Z0' * (Ndiv - diag(diag(Ndiv))) * Z0' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%