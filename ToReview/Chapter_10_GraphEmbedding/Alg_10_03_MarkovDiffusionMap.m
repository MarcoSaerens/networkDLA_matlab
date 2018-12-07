function X = Alg_10_03_MarkovDiffusionMap(A, t, p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Gilissen & Edouard Lardinois, revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the diffusion map embedding of a graph.
%
% INPUT:
% -------
% - A : a (n x n) weighted adjacency matrix associated with an undirected
%       graph G.
% - t:  a integer time step parameter with t >= 1.
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
[n, m] = size(A);
if n ~= m
    error('The adjacency matrix must be a square matrix.')
end

% Check if symmetric matrix / graph is undirected
if ~isequal(A, A')
    error('The adjacency matrix is not symmetric.')
end

% Check if time parameter is greater or equal to 1 
if t < 1
    error('The time step parameter must be greater or equal to 1.')
end

% Check if number of dimensions is pertinent
if (p + 1) > n
    error('The number of kept dimensions must be less than the number of dimensions.');
end

%% Algorithm 

% Compute the diagonal degree matrix
Diag_d = diag(A * ones(n, 1));

% Compute the transition matrix
P = (Diag_d)^(-1) * A;

% Compute the p dominant eigenvectors of K, K * u_k = lambda_k * u_k
% corresponding to positive eigenvalues
[U, Lambda] = eigs(P, p+1);

% Sort the eigenvectors and eigenvalues in decreasing order of eigenvalue
% modulus
[lambda, index] = sort(diag(abs(Lambda)), 'descend');
Lambda = diag(lambda(2:p+1));
U  = U(:, index(2:p+1));

% Stack the coordinate vectors in X
X = U * abs(Lambda^t);

end