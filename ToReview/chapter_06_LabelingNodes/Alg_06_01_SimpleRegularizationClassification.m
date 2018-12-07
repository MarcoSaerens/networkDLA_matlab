function l_hat = Alg_06_01_SimpleRegularizationClassification(A, Y, lambda)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: A simple regularization framework for labeling nodes.
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix, representing a undirected 
%       graph.
% - Y : the (n x m) binary matrix containing label indicator vectors on
%       its columns for m classes.
% - lambda : a strictly positive regularization parameter. 
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

% Check if Y has the right number of rows
[n_Y, m] = size(Y);
if n_Y ~= n
    error('The binary matrix for classes does not correspond to the adjacency matrix.')
end

% Check if lambda is stricly positive
if lambda <= 0
    error('The regularization parameter is not stricly positive');
end

%% Algorithm

% Vectors of ones
e_n = ones(n, 1);
e_m = ones(m, 1);

% A diagonal matrix Gamma inficating which nodes are labeled 
Gamma = diag(Y * e_m);

% The degree matrix
Diag_d = diag(A * e_n);

% The Laplacian matrix
L = Diag_d - A;

% Compute the sum-of-similarities scores for each class
Y_hat_star = linsolve( Gamma + lambda*L, Gamma * Y);

% The class label vector
[~, l_hat]  = max(Y_hat_star, [], 2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
