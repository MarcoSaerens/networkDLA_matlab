function D = Alg_03_05_BagOfHittingPathsFreeEnergyDistance(A, C, theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Maxime Duyck revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the bag-of-hitting-paths potential, or free energy,
%              distance matrix of a weighted and directed graph 
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
% - D : the (n x n) bag-of-hitting-paths free-energy distance matrix 
%       containing the pairwise distances between node
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

% Take elementwise logarithm for computing the potentials
Phi = -log(Z_h) / theta;

% Symmetrize the matrix 
D = (Phi + Phi') / 2;

% Put diagonal to zero
D = D - diag(diag(D));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%