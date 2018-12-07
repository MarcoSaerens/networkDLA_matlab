function l_hat = Alg_06_06_RidgeWithLaplacianRegClassification(A, X, Y, lambda, mu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: A ridge regression with Laplacian regularization for 
%              labeling nodes of a graph and integrating features 
%              available on the nodes.
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix, representing an undirected 
%       graph.
% - X : the (n x q) data matrix containing features vectors on the nodes 
%       on its rows.
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

% Check if X has the right number of rows
[n_X, q] = size(X);
if n_X ~= n
    error('The data matrix does not correspond to the adjacency matrix.')
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

% Vectors of ones, 
e_n = ones(n, 1);
e_m = ones(m, 1);

% A diagonal matrix Gamma inficating which nodes are labeled 
Gamma = diag(Y * e_m);

% The augmented data matrix with a column full of 1's in order to 
% integrate a bias term into the model.
X = [X e_n];

% Identity matrix for size q + 1
I_q = eye(q + 1);

% The row-normalization, or outdegree, matrix. 
Diag_d = diag(A * e_n);

% The Laplacian matrix.
L = Diag_d - A;

% The loop on classes 
Y_hat = zeros(n, m);
for c = 1:m
    % Fit an initial logistic regression 
    w_c_hat = mnrfit(n, 1);
end

% Compute the parameter matrix W_hat for all classes 
W_hat = linsolve( X'*Gamma*X + lambda*X'*L*X + mu*I_q, X'*Gamma*Y );
% Compute the predicted scores for all classes
Y_hat = X * W_hat

% The class label vector
[~, l_hat]  = max(Y_hat, [], 2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
