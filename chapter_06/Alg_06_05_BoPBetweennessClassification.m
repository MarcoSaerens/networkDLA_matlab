function l_hat = Alg_06_05_BoPBetweennessClassification(A, C, Y, theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: The Bag-of-paths group betweennness approach for labeling 
%              nodes.
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix, representing a directed, 
%       strongly connected, aperiodic graph, without any self-loop.
% - C : the (n x n) cost matrix, associated with G.
% - Y : the (n x m) binary matrix containing label indicator vectors on
%       its columns for m classes.
% - theta : the stricly positive inverse temperature parameter.
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

% Check if the cost matrix correspond to A
[n_C, m_C] = size(A);
if (n_C ~= n) || (m_C ~= m)
    error('The cost matrix does not correspond to the adjacency matrix.')
end

% Check if Y has the right number of rows
[n_Y, m] = size(Y);
if n_Y ~= n
    error('The binary matrix for classes does not correspond to the adjacency matrix.')
end

% Check if theta is stricly positive.
if theta <= 0
    error('The inverse temperature parameter is not stricly positive.')
end

%% Algorithm

% Vectors of ones, identity matrix
e = ones(n, 1);
I = eye(n);

% The row-normalization, or outdegree, matrix. 
Diag_d = diag(A * e);

% The transition matrix.
P_ref = Diag_d^(-1) * A;

% The elementwise expenential and multiplication
W = P_ref .* exp(-theta * C);

% The fundamental matrix
Z = (I - W)^(-1);

% Set the diagonal to zero
Z_0 = Z;
Z_0(1:n+1:end) = 0;

% Inverse degree diagonal matrix
Diag_Z_inv = diag(diag(Z))^(-1);

% The loop on m classes
G = zeros(n, m);
for c = 1:m
    y_c = Y(:, c);
    g_c = Diag_Z_inv * ( ((Z_0' * y_c).*(Z_0 * y_c)) - (Z_0' .* Z_0)*y_c );
    G(:, c) = g_c / norm(g_c, 1);
end

% The class label vector
[~, l_hat]  = max(G, [], 2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
