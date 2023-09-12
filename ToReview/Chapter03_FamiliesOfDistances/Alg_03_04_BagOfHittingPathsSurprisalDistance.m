function D_h = Alg_03_04_BagOfHittingPathsSurprisalDistance(Pi_h)

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
% - Pi_h : The (n x n) hitting paths probability matrix, associated to a
%          weighted and directed graph, containing the bag-of-hitting-paths 
%          probabilities (included zero-length paths). To compute it, use 
%          the algorithm 3.3.
% 
% OUTPUT:
% -------
% - D_h : The (n x n) bag-of-hitting-paths surprisal distance matrix
%         containing the pairwise distances between nodes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Checks of arguments 
 
% Check if squared matrix 
[n, m] = size(Pi_h);
if n ~= m 
    error('The probability matrix is not squared')
end

%% Algorithm 

% Take elementwise logarithm 
D_h = -log(Pi_h);

% Symmetrize the matrix 
D_h = (D_h + D_h') / 2;

% Put diagonal to zero
D_h = D_h - diag(diag(D_h));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%