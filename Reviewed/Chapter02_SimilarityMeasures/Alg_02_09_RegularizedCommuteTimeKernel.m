function output = Alg_02_09_RegularizedCommuteTimeKernel(A, alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Maxime Duyck revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the regularized commute-time kernel and the 
%              random-walk with restart similarity.
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix of a undirected graph G.
% - alpha : the regularized parameter between 0 and 1.
%
% OUTPUT:
% -------
% - output : a structure containing :
%           output.K_RCT  - the (n x n) regularized commute-time kernel matrix.
%           output.K_RWR - the (n x n) random walk with restart similarity matrix.
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
if (alpha < 0) || (alpha > 1)
    error('The regularized parameter is outside the valid range (should be in [0, 1]).')
end

%% Algorithm 

% Degree matrix 
Diag_d = diag( A*ones(n, 1) );

% Regularized commute-time kernel matrix
output.K_RCT = (Diag_d - alpha*A)^(-1);

% Random walk with restart similarity matrix
output.K_RWR = output.K_RCT * Diag_d;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%