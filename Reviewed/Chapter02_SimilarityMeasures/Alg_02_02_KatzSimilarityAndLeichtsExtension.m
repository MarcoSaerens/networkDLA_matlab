function K = Alg_02_02_KatzSimilarityAndLeichtsExtension(A, alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Maxime Duyck, revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computing the Katz similarity matrix and Leicht's 
%              extension of it, for an undirected weighted graph.
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix of an undirected graph G.
% - alpha : the discounting factor ( 0 < alpha < 1/spectral_radius(A) )
% 
% OUTPUT:
% -------
% - K : a structure containing :
%           K.Katz - the (n x n) Katz similarity matrix.
%           K.Leicht - the (n x n) Leicht's extension - the 
%                      degree-weighted Katz similarity matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments 

% Check if squared matrix 
[n, m] = size(A);
if n ~= m 
    error('The adjacency matrix is not squared.')
end

% Check if symmetric matrix / graph is undirected
if ~isequal(A, A')
    error('The adjacency matrix is not symmetric.')
end

% Check if alpha has a correct value 
inv_spectral_radius = 1/max( eig(A) );
if (alpha < 0) || (alpha > inv_spectral_radius)
    error('The discounting factor is outside the valid range (should be in (0, %d) ).', inv_spectral_radius)
end

%% Algorithm 

% Identity matrix 
I = eye(n);

% Auxiliary matrix
Aux = (I - alpha*A)^(-1); 

% Katz similarity matrix 
K.Katz = Aux - I;

% the invers diagonal degree matrix
Diag_d_inv = diag( A*ones(n,1) )^(-1);

% Leicht's extension 
K.Leicht = Diag_d_inv * Aux * Diag_d_inv;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
