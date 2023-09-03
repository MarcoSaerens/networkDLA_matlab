function K_MRL = Alg_02_07_ModifiedRegularizedLaplacianKernel(A, alpha, gamma)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Maxime Duyck revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the modified regularized Laplacian kernel matrix 
%              of a weighted undirected graph.
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix of a undirected graph G.
% - alpha : the discounting parameter, stricly higher than 0.
% - gamma : the parameter controlling importance and relatedness, must be
%           between 0 and 1 (provides the regularized Laplacian kernel 
%           if gamma = 1).
% 
% OUTPUT:
% -------
% - K_MRL : The (n x n) modified regularized Laplacian kernel matrix.
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

% Check if alpha has a correct value 
if alpha <= 0 
    error('The discounting factor is outside the valid range (should be > 0).')
end

% Check if gamma has a correct value 
if (gamma < 0) || (gamma > 1)
    error('The parameter controlling importance and relatedness  is outside the valid range (should be in [0, 1]).')
end

%% Algorithm 

% Degree matrix 
Diag_d = diag( A*ones(n, 1) );

% Modified Laplacian Matrix 
L_gamma = gamma*Diag_d - A ;

% The modified regularized Laplacian kernel matrix
K_MRL = (eye(n) + alpha*L_gamma)^(-1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
