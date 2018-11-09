function x = Alg_05_06_PageRank(A, u, alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: François Fouss revised by XXXX (2017).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the vector containing the PageRank with
% personalization scores of all the nodes of a directed graph.
%
% INPUT:
% ------- 
% - A : the (n x n) weighted adjacency matrix of a directed, strongly 
% connected graph G, were dangling nodes were made absorbing.
% - u : the n-dimensional nonnegative personalization (column) normalized
% vector.
% - alpha : the parameter to assure regularity of G, with 0 < alpha < 1. 
% 
% OUTPUT:
% -------
% - x : The (n x 1) PageRank score vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments 
 
% Check if A is a squared matrix 
[n, m] = size(A);
if n ~= m
    error('The adjacency matrix is not squared')
end

% Check if dangling nodes were made absorbing
d  = sum(A,2); % Outdegree vector
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

% A vector of ones
e = ones(n,1);

% The diagonal matrix containing the outdegree of the nodes
D  = diag(A*e);

% The transition probability matrix
P = D^-1 * A;

% Initialization of x by indegrees
x = A' * e / sum(abs(A' * e));

% Defining the stop condition for convergence
prec = 0.00001;

% Initializing the stop variable
stop = Inf;

while (stop > prec) % Test of the convergence
    old_x = x;
    
    % Update the score vector x
    x = alpha * P' * x + (1-alpha) * u;
    % Normalize the score vector x
    x = x / sum(x);

    % Update the stop variable
    stop = sum(abs(x - old_x)) / sum(old_x);
end

