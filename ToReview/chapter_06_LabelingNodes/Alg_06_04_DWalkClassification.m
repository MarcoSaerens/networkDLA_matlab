function l_hat = Alg_06_04_DWalkClassification(A, Y, alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: The D-walk approach for labeling the nodes of a graph 
%              without self-loops.
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix, representing a directed, 
%       strongly connected, aperiodic graph, without any self-loop.
% - Y : the (n x m) binary matrix containing label indicator vectors on
%       its columns for m classes.
% - alpha : the killing rate of the ranom walker, controlling the 
%           probability of disappearing at each time step, 
%           with 0 <= alpha <= 1.
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
[n_Y, m] = size(Y);
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

% The transition matrix defining a killed random-walk.
P = alpha * Diag_d^(-1) * A;

% The loop on m classes
N = zeros(n, m);
for c = 1:m
    % The class vector.
    y_c = Y(:, c);
    % The number of nodes in class.
    n_c = sum(y_c);
    % Set rows corresponding to class c to the zero vector.
    Q_c = diag(e - y_c) * P;
    % Compute the group betweenness score for class c
    N(:, c) = linsolve(I - Q_c', 1/n_c * P'*y_c);
end

% The class label vector
[~, l_hat]  = max(N, [], 2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
