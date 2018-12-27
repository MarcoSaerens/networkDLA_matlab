function u = Alg_08_03_KCore(A, k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Marco Saerens (2017).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016),
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the k-core of a weighted, undirected, strongly
% connected, graph G.
%
% INPUT:
% ------
% - A: the nxn adjacency matrix of an unweighted, undirected, strongly
%   connected, graph G containing n nodes.
% - k: the order of the desired core.
%
% OUTPUT:
% ------- 
% - u: a n × 1 membership vector indicating which node is part of the k-core. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check if A is a square matrix 
[n,m] = size(A);
if n ~= m
    display('Error : matrix A must be square');
    return;
end

%% Algorithm
d = sum(A,2); % initial degree of the nodes
u = ones(1,n); % binary indicator vector tracking active nodes
l = 1:n; % current list of labels of nodes

% Loop until the degree of each node is greater or equal to k
while (size(d) > 0) & (min(d) < k)
    [~,i]  = min(d); % index of node with minimum degree which will be removed
    
    A(i,:)  = []; % remove node i (row i) with minimum degree
    A(:,i)  = []; % remove node i (column i) with minimum degree
    u(l(i)) = 0;  % remove node i from the list of active nodes
    l(i) = []; % remove node i from the list of labels
    
    d = sum(A,2); % recompute degree of nodes
end
d'
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
