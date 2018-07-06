function x = Alg_05_04_CitationInfluence(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: François Fouss revised by XXXX (2017).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the vector containing the PageRank with
% personalization scores of all the nodes of a directed graph.
%
% INPUT:
% ------- 
% - A : the (n x n) weighted adjacency matrix of a directed, strongly 
% connected graph G. 
% 
% OUTPUT:
% -------
% - x : The (n x 1) citation influence score vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments 
 
% Check if A is a squared matrix 
[n, m] = size(A);
if n ~= m
    error('The adjacency matrix is not squared')
end

%% Algorithm

% A vector of ones
e = ones(n,1);

% The diagonal matrix containing the outdegree of the nodes
D  = diag(A*e);

% The transition probability matrix, assuming a strongly connected graph
P = D^-1 * A;

% The dominant eigenvectors of the transpose of P
[V,unused] = eigs(P');

% The normalized dominant eigenvector of the transpose of P
pi = V(:,1) / sum(V(:,1));

% The influence score vector
x = D^-1 * pi;

