function [xh, xa] = Alg_05_07_HITS(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: François Fouss revised by XXXX (2017).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes HITS hub and authority scores for all the nodes
% of a directed graph.
%
% INPUT:
% ------- 
% - A : the (n x n) weighted adjacency matrix of a directed, strongly 
% connected graph G.
% 
% OUTPUT:
% -------
% - xh : The (n x 1) hub score vector.
% - xa : The (n x 1) authority score vector.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments 
 
% Check if A is a squared matrix 
[n, m] = size(A);
if n ~= m
    error('The adjacency matrix is not squared')
end

%% Algorithm

% A vector of ones
e = ones(n,1);

% Initialize xh and xa by a uniform vector
xh = e / sqrt(n);
xa = e / sqrt(n);

% Defining the stop condition for convergence
prec = 0.00001;

% Initializing the stop variable
stop = Inf;

while (stop > prec) % Test of the convergence
    old_xh = xh;
    old_xa = xa;
    
    % Update the hub score vector xh
    xh = A * xa / norm(A * xa);
    % Update the authority score vector xa
    xa = A' * xh / norm(A' * xh);

    % Update the stop variable
    stop = sum(abs(xh - old_xh)) / sum(old_xh) + sum(abs(xa - old_xa)) / sum(old_xa);
end

