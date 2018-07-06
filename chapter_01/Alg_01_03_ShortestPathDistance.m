function D = Alg_01_03_ShortestPathDistance(C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Marco Saerens (2014), revised by Ilkka Kivimaki 
%          & Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the shortest path directed distances between all 
%              pairs of nodes by Floyd-Warshall algorithm, in matrix form.
%
% INPUT:
% -------
% - C : the (n x n) nonnegative cost matrix, representing 
%       a directed weighted graph. C(i,j) == Inf <=> A(i,j) == 0 
% 
% OUTPUT:
% -------
% - D : the (n x n) directed shortest path distance matrix between every 
%        pair of nodes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments 

% Check if squared matrix 
[n, m] = size(C);
if n ~= m
    error('The adjacency matrix is not squared.')
end

%% Algorithm

% Column vector full of 1s
e = ones(n,1); 

% Initialize distances to costs
D = C;

% Set diagonal to zero
D(1:n+1:end) = 0; 

% Iterations    
for t = 1:n
     D = min( D, (D(:, t)*(e') + e*D(t, :)) );
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
