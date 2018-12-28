function core = Alg_08_04_CoreNumber(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Marco Saerens (2017).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016),
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the core number of the nodes of an unweighted,
% undirected, strongly connected, graph G.
%
% INPUT:
% ------
% - A: the n x n adjacency matrix of an unweighted, undirected, strongly
%   connected, graph G containing n nodes. It must contain 0s for missing links.
%
% OUTPUT:
% ------- 
% - core: a n × 1 vector containing the core number of each node. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check if A is a square matrix 
[n,m] = size(A);
if n ~= m
    display('Error : matrix A must be square');
    return;
end

%% Algorithm
myMax = 10e+100;
d = sum(A,2); % initial degree of the nodes
core = d; % initialize core numbers to initial degrees
l = 1:n; % current list of indexes of nodes of G to process

% Loop until each node has been processed
while (size(l) > 0)
    [~,i] = min(d(l)); % index of node with minimum degree among list l, which must be processed next
    i1 = l(i); % recover the right index of the current optimum node
    
    for j = find(A(i1,:)) % loop on neighbors of node i1
        if (core(j) > core(i1))
            core(j) = core(j) - 1; % decrease the core number of neighbors
        end
        d(j) = d(j) - 1; % decrease the degree of neighbors j
    end
    l(i) = []; % remove optimal node i from the list of nodes to process
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
