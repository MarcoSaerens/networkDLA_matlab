function [ gbet ] = Alg_04_09_BagOfPathsGroupBetweenness(A, C, hi, hk, theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Bertrand Lebichot, revised by Marco Saerens (2018).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016),
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: The bag of paths approach for computing a group betweenness
% centrality measure on the nodes of a weighted directed, strongly
% connected graph G without self-loops.
%
% INPUT:
% ------
% - A: the nxn adjacency matrix of a weighted directed, stronly
%   connected graph G containing n nodes.
% - C: the nxn cost matrix C associated to G (if not specified,
%   the costs are the inverse of the affinities, but other 
%   choices are possible).
% - hi: a nx1 vector containing binary membership values to group 1 (class
%   i).
% - hk: a nx1 vector containing binary membership values to group 2 (class
%   k).
% - theta: the (scalar) inverse temperature parameter.
%
% OUTPUT:
% ------- 
% - gbet: the nx 1 bag of hitting paths group betweenness vector gbet  
%   containing the probabilities that a randomly chosen path between
%   group 1 and group 2 visits intermediate nodes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check if A is a square matrix 
[n,m] = size(A);
if n ~= m
    display('Error : matrix A must be square');
    return;
end

%% Algorithm
d = A * ones(n,1); % the outdegree vector
Dinv = diag(1./d); % the inverse outdegree square matrix
Pref = Dinv * A;  % the reference transition probabilities matrix

W = Pref .* exp(-theta*C); % computation of the W matrix

I  = eye(n); % identity matrix (same size as W)
Z  = (I - W)\I; % the fundamental matrix Z = inv(I - W)
Z0 = Z - diag(diag(Z)); % set diagonal to zero

Dzinv = diag(diag(Z).^-1); % the inverse of diagonal matrix Dz

% Computation of the BoP group betweenness
gbet = Dzinv * ( (Z0' * hi) .* (Z0 * hk) );
gbet = gbet / norm(gbet,1); % normalize the result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
