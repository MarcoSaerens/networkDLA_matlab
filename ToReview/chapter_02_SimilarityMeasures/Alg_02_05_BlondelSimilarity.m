function K = Alg_02_05_BlondelSimilarity(A, B, threshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Maxime Duyck, revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the Blondel and al. similarity matrix 
%              between nodes of two weighted, directed graphs
%
% INPUT:
% -------
% - A : the (n_A x n_A) weighted adjacency matrix of a directed graph G_A.
% - B : the (n_B x n_B) weighted adjacency matrix of a directed graph G_B.
% - (optional) threshold : the convergence threshold (default = 1e-3).
% 
% OUTPUT:
% -------
% - K : The n_A x n_B similarity matrix between nodes of graphs G_A 
%       and G_B.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments 

% Default convergence threshold
if nargin == 2
    threshold = 1e-3;
end

% Check if A is a squared matrix 
[n_A, m_A] = size(A);
if n_A ~= m_A
    error('The adjacency matrix A is not squared.')
end

% Check if B is a squared matrix 
[n_B, m_B] = size(B);
if n_B ~= m_B
    error('The adjacency matrix B is not squared.')
end

%% Algorithm 

% Initialistion K = E
K = ones(n_B, n_A);

% Iterations
convergence = 0;
i = 0;
while ~convergence
    
    K_old = K;
    
    % Matrix of iterated scores containingsimilarities between nodes
    K = (B*K*A') + (B'*K*A); 
    
    % Normalize the matrix by its Frobenius Norm
    K = K / norm(K,'fro'); 
    
    i = i + 1; 
    if (~mod(i,2)) & (min(abs(K-K_old)) < threshold) % Convergence condition
        convergence = 1; 
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%