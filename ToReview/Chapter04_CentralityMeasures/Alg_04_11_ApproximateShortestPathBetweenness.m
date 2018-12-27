function bet = Alg_04_11_ApproximateShortestPathBetweenness(A, C, theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Ilkka Kivimaki (2015), revised by Marco Saerens (2018).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the RSP betweenness without weighting by the reference
% random walk so that we recover the shortest path betweenness when theta
% becomes large.
%
% INPUT:
% ------ 
% - Pref: the (n x n) adjacency matrix of a weighted, strongly connected,
%   possibly directed graph G.
% - C: the (n x n) cost matrix C associated to G (default = 1./A).
% - theta: the (scalar) non-negative inverse temperature parameter.
% 
% OUTPUT:
% -------
% - bet: the nx 1 randomized shortest paths betweenness vector containing  
%   the frequency that a sampled path visits an intermediate node.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Some preprocessing and checks
[n,m] = size(A);
if n ~= m
    error('Matrix A must be square')
end

eps = 10^(-40); % tolerance level to be set by user

% % Eventually normalize C
% nEdges = sum(sum(A>eps)); % compute number of edges
% c = sum(sum(C(C<(1/eps))))/nEdges; % compute mean cost of edges
% C = C/c; % normalise the cost matrix so that mean cost of edges is 1

%% Algorithm
W = exp(-theta*C); % RSP's W without Pref for recovering shortest-path betweenness when theta is large
I = eye(n); % identity matrix (same size as W)

% If spectral radius is greater or equal to 1, series does not converge
rs = max(sum(W,2)); % check first max row sum
if rs > (1 - eps)
    r = abs(eigs(W,1)); % then compute spectral radius of W
    if r > (1 - eps)
        error('Error: series could not converge: %0.5g',r')
    end
end

Z = (I - W)\I; % the fundamental matrix Z = inv(I - W)

if any(any(Z<=0 | isinf(Z)))
    error('Z contains zero or Inf values - either graph is not strongly connected or beta is too small/large');
end

Zdiv = 1./Z; % the matrix containing elements 1/Z(i,j)

DZdiv = diag(diag(Zdiv)); % diagonal matrix containing the diagonal
                          % elements 1/Z(i,i)

bet = diag( Z * transpose(Zdiv - (n-1)*DZdiv) * Z ) - n*diag(Z); % the RSP betweenness vector
bet = bet/((n-1)*(n-2)); % shortest-path betweenness is usually not normalized so this line must 
                         % be commented if you want to compare with shortest-path betweenness
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%