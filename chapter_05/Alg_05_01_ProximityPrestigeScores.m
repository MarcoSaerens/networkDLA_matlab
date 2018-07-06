function x = Alg_05_01_ProximityPrestigeScores(D)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: François Fouss revised by XXXX (2017).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the vector containing the proximity prestige scores
% of all the nodes of a directed graph.
%
% INPUT:
% ------- 
% - D : the (n x n) directed shortest path distance matrix between every 
%        pair of nodes of a directed weighted graph, computed by applying
%        Alg_01_03_ShortestPathDistance on C, the (n x n) nonnegative cost
%        matrix, representing a directed weighted graph.
% 
% OUTPUT:
% -------
% - x: The (n x 1) proximity prestige scores vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments 
 
% XXXX Je suppose qu'il n'y a pas de check puisqu'on a utilisé l'Alg_01_03

%% Algorithm
% The number of nodes of the graph
[n, unused] = size(D);

% Initialization of the number of nodes that can reach a node
nb = zeros(n,1);

% Initialization of the sum-of-the-distances vector
sum = zeros(n,1);

% Initialization of the proximity prestige score vector
x = zeros(n,1);

for j = 1:n
    for i = 1:n
        if D(i,j) > 0
            delta1 = 1;
        else
            delta1 = 0;
        end
        if D(i,j) < Inf
            delta2 = 1;
            sum(j) = sum(j) + delta2 * D(i,j);
        else
            delta2 = 0;
        end
        nb(j) = nb(j) + delta1 * delta2;
        
    end
    x(j) = (nb(j) / (n - 1)) / (sum(j) / nb(j));
end