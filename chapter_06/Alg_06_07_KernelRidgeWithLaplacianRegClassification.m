function l_hat = Alg_06_07_KernelRidgeWithLaplacianRegClassification(A, K, Y, lambda, mu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: A kernel ridge regression with Laplacian regularization for 
%              labeling nodes of a graph and integrating features 
%              available on the nodes.
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix, representing an undirected 
%       graph.
% - K : the (n x n) kernel matrix containing features-based similarities
%       between nodes.
% - Y : the (n x m) binary matrix containing label indicator vectors on
%       its columns for m classes.
% - lambda : a positive regularization parameter.
% - mu : a positive regularization parameter.
%
% OUTPUT:
% -------
% - l_hat : a (n x 1) class label vector containing the predicted class.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments 

% Check if squared matrix 
[n, m] = size(A);
if n ~= m
    error('The adjacency matrix is not squared.')
end

% Check if symmetric matrix / graph is undirected
if ~isequal(A, A')
    error('The adjacency matrix is not symmetric.')
end

% Check if K corresponds to A
[n_K, m_K] = size(K);
if (n_K ~= n) || (m_K ~= m)
    error('The kernel matrix does not correspond to the adjacency matrix.')
end

% Check if Y has the right number of rows
[n_Y, m] = size(Y);
if n_Y ~= n
    error('The binary matrix for classes does not correspond to the adjacency matrix.')
end

% Check if regulariztation parameters are positive.
if (lambda < 0) || (mu < 0) 
    error('One of the regularization parmaters is not positive.')
end

%% Algorithm

% Vectors of ones, identity matrix
e_n = ones(n, 1);
e_m = ones(m, 1);
I_n = eye(n);

% A diagonal matrix Gamma inficating which nodes are labeled 
Gamma = diag(Y * e_m);

% The row-normalization, or outdegree, matrix. 
Diag_d = diag(A * e_n);

% The Laplacian matrix.
L = Diag_d - A;

% Compute the parameter matrix Beta_hat for all classes 
Beta_hat = linsolve(Gamma*K + lambda*L*K + mu*I_n, Gamma*Y);
% Compute the predicted scores for all classes
Y_hat = K * Beta_hat

% The class label vector
[~, l_hat]  = max(Y_hat, [], 2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
