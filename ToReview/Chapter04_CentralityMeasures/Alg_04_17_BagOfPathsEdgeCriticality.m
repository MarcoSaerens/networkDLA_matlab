function Cr = Alg_04_17_BagOfPathsEdgeCriticality(A, C, theta, Edg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Marco Saerens (2016).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016),
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%         See also: Lebichot B. and Saerens M. (2018), "A bag-of-paths
%         network criticality measure". Neurocomputing, 275, pp. 224-236.
%
% Description: The bag of paths approach for computing an approximate
% criticality measure on the edges of a weighted directed or undirected,
% strongly connected, graph G without self-loops. See the Neurocomputing
% paper for more information.
%
% INPUT:
% ------
% - A: the n x n adjacency matrix of a weighted undirected or directed,
%   strongly connected, graph G containing n nodes.
% - C: the n x n cost matrix C associated to G (usually,
%   the costs are the inverse of the affinities, but other 
%   choices are possible).
% - theta: the (scalar) inverse temperature parameter.
% - Edg: a q x 2 matrix containing the indices of q edges (i,j) as rows
%
% OUTPUT:
% ------- 
% - cr: the q x 1 bag of paths criticality vector containing edge
%   criticalities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check if A is a square matrix 
[n,m] = size(A);
if n ~= m
    display('Error : matrix A must be square');
    return;
end

q = length(Edg);

%% Algorithm

% Initialization
e  = ones(n,1); % column vector full of 1s
I  = eye(n,n); % identity matrix (same size as A)
Cr = zeros(n,n); % initialize edge criticality matrix
eps = 10e-50; % precision

do = A * e; % the outdegree matrix
di = A' * e; % the indegree matrix
Dinv = diag(1./do); % inverse of diagonal outdegree matrix 

W = exp(-theta*C); % the auxiliary matrix W in the book
% Pref = Dinv * A; % the reference transition probabilities matrix
% W = Pref .* exp(-theta*C); % an alternative: the auxiliary matrix W in the Neurosomputing paper

sr = abs( eigs(W,1,'LM') ); % compute spectral radius of matrix W
if sr < (1 - eps)
    Z = (I-W)\I; % compute the fundamental matrix
    for row = 1:length(Edg) % loop on edges; compute criticality on each edge in turn
        i = Edg(row,1); % the currently processed starting node of the edge
        j = Edg(row,2); % the currently processed ending node of the edge
        if (do(i) > 1) && (di(j) > 1)
            Cr(i,j) = -(1/theta) * real( log( 1 + W(i,j) * Z(j,i) - (Z(i,i)*W(i,j)*Z(j,j))/(Z(i,j)) + realmin ) );
        end
    end
else
    error('theta is too small so that the spectral radius of W is larger than 1');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
