function x = Alg_05_06_PageRank(A, u, alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: François Fouss revised by Marco Saerens (2017).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the vector containing the PageRank with
% personalization scores of all the nodes on a directed, weighted,
% graph.
%
% INPUT:
% ------ 
% - A: the (n x n) adjacency matrix of a weighted directed graph G,
% where dangling nodes were made absorbing.
% - u: the n-dimensional nonnegative personalization (column) vector,
% normalized to sum to one.
% - alpha: the parameter to assure regularity of G, with 0 < alpha < 1. 
% 
% OUTPUT:
% -------
% - x: the (n x 1) PageRank score vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Check if A is a square matrix 
[n, m] = size(A);
if n ~= m
    error('The adjacency matrix is not square')
end

% Check if dangling nodes were made absorbing
d  = sum(A,2); % outdegree vector
if min(d) == 0
    error('The dangling nodes were not all made absorbing')
end

% Check if u is normalized 
if abs(sum(u) - 1) > 0.00001
    error('The personalization vector is not normalized')
end

% Check if alpha is strictly positive and less than 1
if alpha <= 0 || alpha >= 1 
    error('The parameter to assure regularity of G must be between 0 and 1 (excluded)')
end


%% Algorithm
e = ones(n,1); % a vector of ones

% The diagonal matrix containing the outdegree of the nodes
D = diag(A*e);

P = D^-1 * A; % the transition probability matrix
x = A' * e / sum(abs(A' * e)); % initialization of x by normalized indegrees

prec = 0.000001; % defining the stop condition for convergence
stop = Inf; % initializing the stop variable monitoring convergence

while (stop > prec) % test of convergence
    old_x = x;
    
    % Update the score vector x
    x = alpha * P' * x + (1-alpha) * u;
    % Normalize the score vector x (not really needed; included for safety)
    x = x / sum(x);

    % Update the stop variable
    stop = sum(abs(x - old_x)) / sum(old_x);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
