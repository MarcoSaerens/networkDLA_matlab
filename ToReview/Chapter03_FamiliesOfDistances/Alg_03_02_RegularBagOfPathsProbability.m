function Pi = Alg_03_02_RegularBagOfPathsProbability(A, C, theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Maxime Duyck revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the regular bag-of-paths probability matrix 
%              (including zero-length paths).
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
% - Pi: The (n x n) bag-of-paths probability matrix containing the
%       probability of drawing a path starting in node i and ending 
%       in node j from a bag of paths, when sampling paths according to 
%       a Gibbs-Boltzmann distribution.
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

% The vector of ones and the identity matrix
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

% The normalization factor - the partition function 
z_dotdot = e'*Z*e;

% The regular bag-of-paths probability matrix
Pi = Z / z_dotdot;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
