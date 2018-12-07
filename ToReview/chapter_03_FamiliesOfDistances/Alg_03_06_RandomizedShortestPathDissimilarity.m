function D = Alg_03_06_RandomizedShortestPathDissimilarity(A, C, theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Maxime Duyck revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the randomized shortest-path dissimilarity matrix  
%              for hitting paths.
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix of an directed, strongly 
%       connected graph G. 
% - C : the (n x n) cost matrix associated to G. 
% - theta : the stricly positive inverse temperature parameter. 
% 
% OUTPUT:
% -------
% - D : the (n x n) randomized shortest-path dissimilarity matrix.
%       containing the pairwise distances between nodes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments 
 
% Check if A is a squared matrix 
[n, m] = size(A);
if n ~= m
    error('The adjacency matrix is not squared')
end

% Check if C has the same dimension as A 
[o, p] = size(C);
if (o ~= n) || (p ~= m)
    error('The cost matrix does not correspond to the adjacency matrix')
end

% Check if theta is positive
if theta <= 0
    error('The inverse temperature parameter is not stricly positive')
end

%% Algorithm 

% The vector of ones and the idendity matrix
e = ones(n, 1);
I = eye(n); 

% The degree matrix
Diag_d = diag(A*e);

% The reference transition probability matrix 
P_ref = Diag_d^(-1) * A;

% The killed Markov chain matrix 
W = P_ref .* exp(-theta*C);

% The fundamental matrix 
Z = (I - W)^(-1);

% The expected costs for nonhitting paths
CdotW = (C .* W);
CdotW(isnan(CdotW)) = 0;
S = (Z*CdotW*Z) ./ Z;

% The randomized shortest path dissimilarity matrix for hitting paths
D = (S + S' - e*(diag(S))' - diag(S)*e') / 2;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
