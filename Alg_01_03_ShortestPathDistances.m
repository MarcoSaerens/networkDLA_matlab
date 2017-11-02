%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = Alg_01_03_ShortestPathDistances(C)
% Authors: Marco Saerens (2014),revised by Ilkka Kivimaki (2017).
% Direct source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
% "Algorithms and models for network data and link analysis".
% Cambridge University Press.
%
% Computes the shortest path directed distances between all pairs of nodes
% by Floyd-Warshall algorithm, in matrix form.
%
% INPUT:
% ------
% - A weighted directed strongly connected graph without self-loop.
% - C is a n x n non-negative cost matrix with cij=infinity for missing links.
% 
% OUTPUT:
% -------
% - D, the n x ndirected shortest path distance matrix between every pair of nodes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,n] = size(C); % square matrix
e = ones(n,1); % column vector full of 1s

D = C; % initialize distances to costs
D(1:n+1:end) = 0; % set diagonal to zero
    
for k=1:n
     D = min( D, (D(:,k)*(e') + e*D(k,:)) )
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
