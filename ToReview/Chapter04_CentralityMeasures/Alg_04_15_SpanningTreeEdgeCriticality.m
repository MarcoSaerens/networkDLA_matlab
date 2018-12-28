function Cr = Alg_04_15_SpanningTreeEdgeCriticality(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Nassim Benoussaid, revised by Marco Saerens (2018).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016),
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Edge Criticality measure based on the spanning-trees
% formalism.
%
% INPUT:
% ------
% - A: the n x n adjacency matrix of a weighted, undirected, strongly
%   connected, graph G containing n nodes.
%
% OUTPUT:
% ------- 
% - Cr: the nx n spanning trees criticality matrix containing edge
%   criticalities, and corresponding to the number of times that the
%   edge lies on a sampled spanning tree.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check if A is a square matrix 
[n,m] = size(A);
if n ~= m
    display('Error : matrix A must be square');
    return;
end

%% Some initializations
e  = ones(n,1); % vector full of 1s
I  = eye(n); % identity matrix (same size as A)
Cr = zeros(n,n); % initialize the criticality matrix

%% Algorithm
L = diag(A'*e) - A; % Laplacian matrix computed from in-degrees (convention used in matrix-tree theorem)

L1 = L;
L1(1,1) = L1(1,1) + 1; % augmented Laplacian matrix
delta = det(L1); % compute the determinant. Not mandatory as it corresponds
% to a simple rescaling factor (total number of spanning trees)

Z  = L1\I; % the inverse of L1 matrix

Cr = delta * (A .* (e*(diag(Z))' - Z')); % compute the directed edge criticalities in matrix form
Cr(:,1) = 0; % for safety, set criticality to 0 on first column, j=1. Not mandatory
% as the first column should automatically be equal to zero (apart from
% numerical approximation errors)

Cr = (Cr + Cr'); % symmetrize the criticality matrix in case of undirected graph

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%