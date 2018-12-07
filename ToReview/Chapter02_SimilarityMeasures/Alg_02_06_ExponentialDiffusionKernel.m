function K = Alg_02_06_ExponentialDiffusionKernel(A, alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Maxime Duyck revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the exponential diffusion kernel and the 
%              Laplacian exponential diffusion kernel matrix of a 
%              connected weighted undirected graph 
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix of a undirected, connected 
%       graph G.
% - alpha : the discounting parameter, stricly higher than 0.
% 
% OUTPUT:
% -------
% - K : a structure containing :
%           K.ED  - the (n x n) exponential diffusion kernel matrix, based
%                   on matrix exponential.
%           K.LED - the (n x n) Laplacian exponential diffusion kernel 
%                   matrix commute-time distances.
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
if alpha <= 0 
    error('The discounting factor is outside the valid range (should be > 0)')
end

%% Algorithm 

% Degree Matrix
Diag_d = diag( A*ones(n,1) );

% Laplacian Matrix 
L = Diag_d - A; 

% Exponential diffusion kernel matrix 
K.ED = expm( alpha*A );

% Laplacian Exponential diffusion kernel matrix 
K.LED = expm( -alpha*L );

end

