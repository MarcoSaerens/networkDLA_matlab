function D = Alg_01_02_ShortestPathDistanceElementwise(C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Marco Saerens (2014), revised by Ilkka Kivimaki
%          & Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis".
%         Cambridge University Press.
%
% Description: Computes the shortest path directed distances between all
%              pairs of nodes by Floyd-Warshall algorithm, in elementwise form.
%
% INPUT:
% -------
% - C : the (n x n) nonnegative cost matrix, representing
%       a directed weighted graph. C(i,j) == Inf <=> A(i,j) == 0.
%
% OUTPUT:
% -------
% - D : the (n x n) directed shortest path distance matrix between every
%       pair of nodes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments

% Check if squared matrix
[n, m] = size(C);
if n ~= m
    error('The adjacency matrix is not squared.')
end

%% Algorithm

% Initialize distances to costs
D = C;

% Set diagonal to zero
D(1:n+1:end) = 0;

% Iterations
for t = 1:n % enumerate intermediate nodes
    for i = 1:n % enumerate starting nodes
        for j = 1:n % enumerate ending nodes
            D(i,j) = min( D(i,j), (D(i,t) + D(t,j)) ); % recompute distances
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
