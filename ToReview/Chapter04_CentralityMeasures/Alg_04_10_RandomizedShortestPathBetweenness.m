function bet = Alg_04_10_RandomizedShortestPathBetweenness(Pref, C, theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Ilkka Kivimaki (2015), revised by Marco Saerens (2018).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the simple RSP betweenness as the expected number
% of visits to each node over RSP distributions based on values of theta.
%
% INPUT:
% ------ 
% - Pref: the (n x n) transition probability matrix of a weighted, strongly
%   connected, possibly directed graph G.
% - C: the (n x n) cost matrix C associated to G (default = 1./A).
% - theta: the (scalar) non-negative inverse temperature parameter.
% 
% OUTPUT:
% -------
% - bet: the RSP betweenness vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Some checks 
[m,n] = size(Pref);
if m ~= n
    error('Matrix Pref must be square')
end

%% Algorithm
W = Pref .* exp(-theta*C); % RSP's W
I = eye(n); % identity matrix (same size as W)

Z = (I - W)\I; % the fundamental matrix Z = inv(I - W)
  
if any(any(Z<=0 | isinf(Z)))
    error('Z contains zero or Inf values - either graph is not strongly connected or beta is too small/large');
end

Zdiv = 1./Z; % the matrix containing elements 1/Z(i,j)

DZdiv = diag(diag(Zdiv)); % diagonal matrix containing the diagonal
                          % elements 1/Z(i,i)

% With this form, RSPBC increase for i = s and i = t and converges to
% the stationary distribution (up to a constant multiplying factor):
bet = (1/((n-1)*n)) * diag( Z * transpose(Zdiv - n*DZdiv) * Z ); % the RSP betweenness vector
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%