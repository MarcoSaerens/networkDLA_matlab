function x = Alg_05_02_Bonacich(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: François Fouss revised by XXXX (2017).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the vector containing the Bonacich scores
% of all the nodes of a directed graph.
%
% INPUT:
% ------- 
% - A : the (n x n) weighted adjacency matrix of a directed, strongly 
%       connected graph G.
% 
% OUTPUT:
% -------
% - x : The (n x 1) vector holding Bonacich's spectral measure of prestige.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments 
 
% Check if A is a squared matrix 
[n, m] = size(A);
if n ~= m
    error('The adjacency matrix is not squared')
end

%% Algorithm

% The dominant eigenvectors of the transpose of A
[V,D] = eigs(A');

% The dominant eigenvector of the transpose of A (i.e., the left
% dominant eigenvector of A, which is defined in absolute value 
x = abs(V(:,1));