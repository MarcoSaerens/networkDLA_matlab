function X = Alg_10_01_KernelPrincipalComponentAnalysis(K, p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Gilissen & Edouard Lardinois, revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes a kernel principal-component analysis of a graph 
%              from a kernel matrix.
%
% INPUT:
% -------
% - K : a (n x n) semidefinite matrix containing the similarities
%       between nodes of graph G.
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
[n, m] =  size(K);
if n ~= m
    error('The Kernel Matrix must be a square matrix.');
end

% Check if number of dimensions is pertinent
if p > n
    error('The number of kept dimensions must be less or equal to the number of dimensions.');
end

%% Algorithm 

% Compute the centering matrix
H = eye(n) - ones(n)/n;

% Center the kernel if not already centered
K = H*K*H; 

% Compute the p dominant eigenvectors of K
[U, Lambda] = eigs(K, p);

% Sort the eignevectors and eigenvalues in decreasing order of eigenvalue
[lamdba, Ind] = sort(diag(Lambda), 'descend');
Lambda = diag(lamdba);
U = U(:, Ind);

% Stack the coordinate vectors in X
X = U * sqrt(Lambda);

end