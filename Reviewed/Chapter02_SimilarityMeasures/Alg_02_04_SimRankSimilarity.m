function K = Alg_02_04_SimRankSimilarity(A, alpha, threshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Maxime Duyck, revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis".
%         Cambridge University Press.
%
% Description: Computes the SimRank similarity matrix between
%              nodes of a directed, weighted, graph
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix of a directed,
%       connected, graph G.
% - alpha : a parameter in [0, 1].
% - (optional) threshold : the convergence threshold (default = 1e-5).
%
% OUTPUT:
% -------
% - K : the (n x n) SimRank similarity matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments

% Default convergence threshold
if nargin == 2
    threshold = 1e-5;
end

% Check if squared matrix
[n, m] = size(A);
if n ~= m
    error('The adjacency matrix is not squared.')
end

% Check if alpha has a correct value
if (alpha < 0) || (alpha > 1)
    error('The discounting factor is outside the valid range (should be in [0, 1]).')
end

%% Algorithm

% Initialistion K = I
K = eye(n);

% Indegree vector
d = A'*ones(n,1);

% Initialization of Q
Q = zeros(n);
for j = 1:n
    if d(j) > 0
        Q(:,j) = A(:,j) / d(j);
    end
end

% Iterations
convergence = 0;
while ~convergence
    K_old = K;
    K = alpha*Q'*K*Q; % update the similarity matrix
    K(1:n+1:end) = (d > 0); % set diagonal to 1 if node has a predecessor
    if max(max( abs(K - K_old) )) < threshold % convergence condition
        convergence = 1;
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%