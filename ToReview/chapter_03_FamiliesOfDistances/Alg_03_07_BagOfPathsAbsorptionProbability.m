function B = Alg_03_07_BagOfPathsAbsorptionProbability(A, C, theta, abs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Maxime Duyck revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the bag-of-paths absorption probabilities to a set
%              of absorbing nodes of a weighted and directed graph  
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix of an directed, strongly 
%       connected graph G. 
% - C : the (n x n) cost matrix associated to G. 
% - theta : the stricly positive inverse temperature parameter.
% - abs : the (n x 1) set of absorbing nodes vector,
%         abs(i) == 1 if i is an absorbing node,
%         abs(i) == 0 otherwise.  
% 
% OUTPUT:
% -------
% - B : the ( (n - |abs|) x |abs| ) matrix containing the bag-of-paths 
%       absoption probabilities for each transient node.
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

% Check if abs has the correct length
if length(abs) ~= n
    error('The set of absorbing nodes do not contain all nodes')
end

% Check if binary vector 
if ~isequal(abs, (abs > 0))
    error('The set of absorbing nodes is not binary')
end

%% Algorithm 

% Initialisation 
e = ones(n,1);
I = eye(n); 

% The set of transient nodes (trans = V\abs)
trans = ones(n,1) - abs;

% The row normalization matrix 
Diag_d = diag(A*e);

% The reference transition probability matrix 
P_ref = Diag_d^(-1) * A;

% Elemental exponential and multiplication 
W_a = P_ref .* exp(-theta*C);

% Set rows corresponding to absorbing nodes to zero
W_a = W_a .* trans;

% The fundamental matrix 
Z_a = (I - W_a)^(-1);

% Extract sub-matrix whose rows correspond to transient nodes and columns
% to absorbing nodes 
B = Z_a(logical(trans), logical(abs));

% Compute normalization factor for each row
Diag_d_B = diag( B * ones(sum(abs),1) );
 
% Normalize the matrix for computing absoption probabilities 
B = Diag_d_B^(-1) * B;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%