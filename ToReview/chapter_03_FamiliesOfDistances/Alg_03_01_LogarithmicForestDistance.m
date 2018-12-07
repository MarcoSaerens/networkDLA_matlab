function D = Alg_03_01_LogarithmicForestDistance(A, alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Authors: Maxime Duyck revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the logarithmic forest distance matrix between 
%              nodes of an undirected graph
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix of an undirected graph G. 
% - alpha: a stricly positive parameter.
% 
% OUTPUT:
% -------
% - D: The (n x n) logarithmic forest distance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments 
 
% Check if squared matrix 
[n, m] = size(A);
if n ~= m
    error('The adjacency matrix is not squared')
end

% Check if symmetric matrix / graph is undirected
if ~isequal(A, A')
    error('The adjacency matrix is not symmetric.')
end

% Check if alpha has a correct value 
if alpha <= 0 
    error('The discounting factor is not stricly positive')
end

%% Algorithm 

% The vector of ones and the idendity matrix
e = ones(n,1);
I = eye(n); 

% The degree matrix
Diag_d = diag(A*e);

% The Laplacian matrix 
L = Diag_d - A;

% Regularized Laplacian Kernel Matrix 
K_RL = (I + alpha*L)^(-1);

% Logarithmic transformation of K_RL 
if alpha == 1
    S = log(K_RL);
else
    S = (alpha - 1)*(log(K_RL)/log(alpha));
end

% Logarithmic forest distance matrix 
D = diag(S)*e' + e*(diag(S))' - 2*S;

end
