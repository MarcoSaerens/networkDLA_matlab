function l_hat = Alg_06_02_HarmonicFunctionClassification(A, Y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: The Harmonic function approach for labeling nodes.
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix, representing a undirected 
%       graph. The node are sorted in such a way that the l first nodes
%       are labeled while the (n - l) last nodes are unlabeled. 
% - Y : the (n x m) binary matrix containing label indicator vectors on
%       its columns for m classes.
%
% OUTPUT:
% -------
% - l_hat : a ( (n - l) x 1) class label vector containing the predicted 
%           class for each unlabeled node.
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
[n_Y, ~] = size(Y);
if n_Y ~= n
    error('The binary matrix for classes does not correspond to the adjacency matrix.')
end

% Check if only the first l nodes are labelled
colSums_Y = sum(Y, 2);
if ~isequal(sort(colSums_Y, 'descend'), colSums_Y)
    error('Only the first l nodes must be labelled (or there are nodes in multiple classes).')
end

%% Algorithm

% Number of labeled nodes
l = sum(colSums_Y);

% Vectors of ones
e = ones(n, 1);

% The degree matrix
Diag_d = diag(A * e);

% The Laplacian matrix
L = Diag_d - A;

% Extract the rows corresponding to unlabeled nodes and the columns 
% corresponding to labeled nodes from the Laplacian matrix
L_ul = L((l + 1):n, 1:l);

% Extract the rows and columns corresponding to unlabeled nodes
L_uu = L((l + 1):n, (l + 1):n);

% Extract the rows corresponding to labeled nodes in Y
Y_l = Y(1:l, :);

% Compute the scores for each class
Y_hat_star = linsolve( L_uu, -L_ul*Y_l);

% The class label vector
[~, l_hat]  = max(Y_hat_star, [], 2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
