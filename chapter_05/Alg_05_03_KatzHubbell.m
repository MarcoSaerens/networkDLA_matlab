function [xk, xh] = Alg_05_03_KatzHubbell(A, u, alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: François Fouss revised by XXXX (2017).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes Katz and Hubbell scores for all the nodes of a
% directed graph.
%
% INPUT:
% ------- 
% - A : the (n x n) weighted adjacency matrix of a directed, strongly 
% connected graph G.
% - u : the n-dimensional nonnegative personalization (column), supposed to
% be normalized, vector (for Hubbell score only).
% - alpha : a discounting parameter, with 0 < alpha < (1 / the spectral
% radius of A). 
% 
% OUTPUT:
% -------
% - xk : The (n x 1) vector holding Katz scores.
% - xh : The (n x 1) vector holding Hubbell scores.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments 
 
% Check if A is a squared matrix 
[n, m] = size(A);
if n ~= m
    error('The adjacency matrix is not squared')
end

% Check if u is normalized 
if abs(sum(u) - 1) > 0.00001
    error('The personalization vector is not normalized')
end

% The largest eigenvalues of A
D = eigs(A);
% Check if alpha is positive and less than the spectral radius of A
if alpha <= 0 || alpha >= 1/abs(D(1)) 
    error('The discounting parameter must be between 0 and the spectral radius of A (excluded)')
end

%% Algorithm

% A vector of ones
e = ones(n,1);

% The identity matrix
I = eye (n);

% The Katz score vector
xk = (alpha * A * (I - alpha * A) ^ -1) * e;

% The Hubbell score vector
xh = ((I - alpha * A) ^ -1) * u;

