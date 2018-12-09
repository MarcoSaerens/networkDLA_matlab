function bet = Alg_04_02_Brandes(C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Masashi Shimbo (2018).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the vector containing the Freeman's centrality scores
%              using Brandes's algorithm.
%
% INPUT:
% ------- 
% - C: the (n x n) cost matrix of a directed graph.
%      To allow a sparse matrix to represent C, components 0 in C are taken
%      as Inf; that is, corresponding edges are nonexistent.
%
% OUTPUT:
% -------
% - bet: The (n x 1) vector of betweenness scores.
%
% NOTE: Requires PriorityQueue.m.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n, nc] = size(C);
if (n ~= nc)
    error('Cost matrix not square');
end

% Both 0 and +Inf in cost matrix C represent "no edge".
% C(C == Inf) = 0; % use 0 uniformly to represent "no edge" from now on.
  
bet = zeros(n, 1);
for i = 1:n
    % pass 1: run Dijkstra's method with node i as the origin
    [~, sigma, pred, S] = Dijkstra(C, i);
    % pass 2: collect dependency scores
    bet = collectDeps(bet, S, sigma, pred, i, n);
end

% N.B. According to the original definition of Freeman's betweenness,
% 'bet' must be halved for undirected graphs, but this process is omitted.

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dist, sigma, pred, S] = Dijkstra(C, origin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: A veriation of Dijkstra's shortest-path algorithm that
%              computes all shortest paths from a given origin node.
% 
% INPUT:
% ------- 
% - C : the (n x n) cost matrix of a directed graph
%     : A 0 entry represents no edge exists between the corresponding nodes
%     : (equivalent of +Inf)
% - origin : the initial node
% 
% OUTPUT:
% -------
% - dist : (1 x n) vector of shortest distance to each node from origin
% - sigma: (1 x n) vector holding the number of shortest-distance paths
%          to each node
% - pred : (n x n) binary matrix representing the shortest-path DAG
%          pred(i, j) == 1 implies that node j is a predecessor of
%          node i on a shortest path from the origin
% - S    : list of nodes reachable from the origin, in decscending order
%          of the distance from the origin (i.e., furthest first, origin
%          last)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = size(C, 1); % number of nodes

% pred(k, :) = vector indicating the optimal predecessors of node k
pred = zeros(n, n); % pred = sparse(n, n); % for huge graphs, use sparse

sigma = zeros(1, n);
sigma(origin) = 1;

% '\delta' on the book
dist = Inf(1, n);
dist(origin) = 0;

S = []; % stack of "closed" nodes

Q = PriorityQueue(n);
Q.insert(origin, 0);

while Q.size ~= 0
    [j, ~] = Q.extractMin;
    S = [j, S];
    % C(:, j) = 0; % optional; not sure how effective this is for speed up
    
    % "relax" edges emanating from j
    succs = (C(j, :) > 0 & C(j, :) ~= Inf);
    altDist = C(j, :) + dist(j) * succs;
    
    updatedOrTied = succs & (altDist <= dist);
    updated = updatedOrTied & (altDist < dist);
    
    for k = find(updated)
        sigma(k) = 0;
        pred(k, :) = 0;
        if dist(k) == Inf
            Q.insert(k, altDist(k));
        else
            Q.decreaseKey(k, altDist(k));
        end
        dist(k) = altDist(k);
    end
  
    pred(updatedOrTied, j) = 1;
    sigma = sigma + sigma(j) * updatedOrTied;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bet = collectDeps(bet, S, sigma, pred, origin, n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Collect the dependency scores of nodes, traversing the
%              shortest-path DAG backward towards the origin.
%
% INPUT:
% ------- 
% - bet   : (n x 1) vector holding the betweenness scores of nodes.
% - sigma : (1 x n) vector holding the number of shortest-distance paths
%           to each node.
% - pred  : (n x n) matrix holding the set of predecessors along shortest
%           paths representing a shortest-path DAG from the origin.
% - origin: the initial node.
% - S     : list of nodes reachable from the origin, in decscending order
%           of the distance from the origin (i.e., furthest node first).
% 
% OUTPUT:
% -------
% - bet: updated vector of betweeness scores (n x 1).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dep = zeros(1, n);
pred(:, origin) = 0;
for k = S
    preds = find(pred(k, :));
    dep(preds) = dep(preds) + (sigma(preds) / sigma(k)) * (1 + dep(k));
    bet(k) = bet(k) + dep(k);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
