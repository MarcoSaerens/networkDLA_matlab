function x = Alg_05_05_LeastSquaresRating(A, Edges)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: François Fouss revised by XXXX (2017).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the vector containing the least squares rating or
% prestige score for all the nodes of a directed graph.
%
% INPUT:
% ------- 
% - A : the (n x n) weighted adjacency matrix of a directed graph G. In the
% language of a competition, aij contains the score of j and aji the score
% of i in the game between i and j.
% - Edges : the set of directed edges of G, assumed reciprocated. Edges is
% represented by a (t x 2) matrix where each line contains a pair of nodes
% connected by a directed edge (source -> destination).
% 
% OUTPUT:
% -------
% - x : The (n x 1) vector holding the ratings or strengths
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments 
 
% Check if A is a squared matrix 
[n, m] = size(A);
if n ~= m
    error('The adjacency matrix is not squared')
end

% Check if Edges contains a source and a destination node columns 
[s, c] = size(Edges);
if c ~= 2
    error('The set of directed edges does not contain 2 columns')
end

%% Algorithm

% A vector of ones
e = ones(n,1);

% Initialize design matrix with zeros
B = zeros(s/2,n);

% Initialize margin vector with zeros
r = zeros((s/2)+1,1);

t = 1;

% Loop on all the reciprocated links
for k = 1:s
    if Edges(k,1) < Edges(k,2)
        % Build the design matrix B
        B(t,Edges(k,1)) = B(t,Edges(k,1)) + 1;
        B(t,Edges(k,2)) = B(t,Edges(k,2)) - 1;
        % Compute the net margin score
        r(t) = A(Edges(k,2),Edges(k,1)) - A(Edges(k,1),Edges(k,2));
        t = t + 1;
    end
end

% Add a row vector full of 1's to B
B = [B;e'];

% Compute the rating/strength vector
x = (B' * B)^-1 * B' * r;


