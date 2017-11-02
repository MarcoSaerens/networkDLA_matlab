function Pi_h = Alg_03_03_BagOfHittingPathsProbability(A, C, theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Maxime Duyck revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the bag-of-hitting-paths probability matrix.
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
% - Pi_h: a structure containing :
%         Pi_h.std = The (n x n) bag-of-hitting-paths probability matrix 
%                    containing the probability of drawing a path starting 
%                    in node i and ending in node j from a bag of paths, 
%                    when sampling paths according to a Gibbs-Boltzmann 
%                    distribution.
%         Pi_h.bar = The (n x n) bag-of-hitting-paths probability matrix 
%                    with zero-length paths.
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

% The column-normalization matrix for the bag-of-hitting-paths 
% probabilities 
Diag_Z = diag(diag(Z));

% Column-normalize the fundamental matrix 
Z_h = Z * (Diag_Z^(-1));

% Compute normalization factor - the partition function 
z_h_dotdot = e'*Z_h*e;

% The bag of hitting paths probability matrix with zero paths included 
Pi_h.std = Z_h/z_h_dotdot;

% Column-normalize the fundamental matrix and substract zero-length
% path contributions
Z_h_bar = Z_h - I;

% Compute normalization factor - the partition function 
z_h_bar_dotdot = e'*Z_h_bar*e;

% The bag of hitting paths probability matrix with zero paths excluded 
Pi_h.bar = Z_h_bar/z_h_bar_dotdot;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

