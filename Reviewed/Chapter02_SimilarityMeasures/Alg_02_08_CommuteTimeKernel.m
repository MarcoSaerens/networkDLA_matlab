function K = Alg_02_08_CommuteTimeKernel(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Maxime Duyck revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the commute-time kernel matrix of a weighted 
%              undirected graph (the usual and the corrected forms).
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix of a undirected graph G.
%
% OUTPUT:
% -------
% - K : a structure containing :
%           K.CT  - the (n x n) commute-time kernel matrix.
%           K.CCT - the (n x n) corrected commute-time kernel matrix.
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

%% Algorithm 

% Utilities
e = ones(n, 1);
E = ones(n);
I = eye(n);
E_div_n = E/n;

% Degree vector, diagonal degree matrix and inverse sqrt degree matrix
d = A*e;
Diag_d = diag(d);
Diag_d_invSqrt = diag(1 ./ sqrt(d));

% Volume of the graph 
vol = sum(d);

% Laplacian Matrix
L = Diag_d - A ;

% Pseudoinverse of the Laplacian matrix
L_plus = (L + E_div_n)^(-1) - E_div_n;

% Commute-time kernel matrix
K.CT = L_plus;

% Centering matrix
H = I - E_div_n;

% "Standardized" modularity matrix
M = Diag_d_invSqrt * ( A - (d*d'/vol) ) * Diag_d_invSqrt; 

% Corrected commute-time kernel matrix
K.CCT = H * Diag_d_invSqrt * M * (I - M)^(-1) * M * Diag_d_invSqrt * H;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%