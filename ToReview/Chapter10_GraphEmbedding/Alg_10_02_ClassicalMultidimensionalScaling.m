function X = Alg_10_02_ClassicalMultidimensionalScaling(D, p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Gilissen & Edouard Lardinois, revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes a classical multidimensional scaling based on 
%              a distance matrix between nodes.
%
% INPUT:
% -------
% - D : a (n x n) symmetric distance matrix between nodes of a graph G.
% - p : a integer containing the number of dimensions kept for 
%       the embedding.
%
% OUTPUT:
% -------
% - X : the (n x p) data matrix containing the coordinates of the nodes 
%       for the embedding.    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments

% Check if squared matrix
[n, m] = size(D);
if n ~= m
    error('The distance matrix must be a square matrix.')
end

% Check if number of dimensions is pertinent
if p > n
    error('The number of kept dimensions must be less or equal to the number of dimensions.');
end

eps = 10e-8; % precision for relevant eigenvalues

%% Algorithm 

% Compute the centering matrix
H = eye(n) - ones(n)/n;

% Elementtwise multiplication - compute squared dissimilarities
D_2 = D .* D;

% Compute the inner products matrix from squared distances
K = -(1/2) * H * D_2 * H;

% Compute the p dominant eigenvectors of K
[U, Lambda] = eigs(K, p);

% Sort the eignevectors and eigenvalues in decreasing order of eigenvalue
[lambda, index] = sort(diag(Lambda), 'descend');
lambda(lambda < eps) = 0; % remove negative eigenvalues - corresponding coordinates will be 0
Lambda = diag(lambda);
U = U(:,index);

% Stack the coordinate vectors in data matrix X (a n x p matrix)
X = real( U * sqrt(Lambda) );

end