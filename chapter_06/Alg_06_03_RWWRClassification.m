function l_hat = Alg_06_03_RWWRClassification(A, Y, alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: The Random-walk with restart approach for labeling nodes.
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix, representing a directed, 
%       strongly connected, aperiodic graph.
% - Y : the (n x m) binary matrix containing label indicator vectors on
%       its columns for m classes.
% - alpha : a parameter verifying 0 <= alpha <= 1. 
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

% Check if Y has the right number of rows
[n_Y, ~] = size(Y);
if n_Y ~= n
    error('The binary matrix for classes does not correspond to the adjacency matrix.')
end

% Check if alpha is in the acceptable range.
if (alpha < 0) || (alpha > 1)
    error('The parameter alpha is not between 0 and 1.')
end

%% Algorithm

% Vectors of ones, identity matrix
e = ones(n, 1);
I = eye(n);

% The row-normalization, or outdegree, matrix. 
Diag_d = diag(A * e);

% The transition matrix
P = Diag_d^(-1) * A;

% Compute the group betweenness scores for each class
Y_hat_star =  linsolve(I - alpha*P, (1 - alpha) * Y * diag(1 ./ sum(Y)) );

% The class label vector
[~, l_hat]  = max(Y_hat_star, [], 2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
