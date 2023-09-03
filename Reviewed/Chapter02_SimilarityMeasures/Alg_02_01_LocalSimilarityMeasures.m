function output = Alg_02_01_LocalSimilarityMeasures(A, i, k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis".
%         Cambridge University Press.
%
% Description: Computes various local similarity indices betweeen
%              nodes i and k.
%
% INPUT:
% -------
% - A : the (n x n) unweighted adjacency matrix of an undirected graph G.
% - i, k : indices of compared nodes.
%
% OUTPUT:
% -------
% - output : a structure containing various scalar similarities:
%              output.sim_direct  - direct similarity.
%              output.sim_common  - common neighbors score.
%              output.sim_pref    - preferential attachment index.
%              output.sim_cos     - cosine coefficient.
%              output.sim_Jaccard - Jaccard index.
%              output.sim_Dice    - Dice coefficient.
%              output.sim_prohub  - hub promoted index.
%              output.sim_dehub   - hub depressed index.
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

% Check if binary matrix
if ~isequal(A, (A > 0))
    error('The adjacency matrix is not binary.')
end

%% Algorithm

p = sum( A(i,:).*A(k,:) ); % common neighbors
q = sum( A(i,:).*(1 - A(k,:)) ); % neighbors only to i
r = sum( (1 - A(i,:)).*A(k,:) ); % neighbors only to k

output.sim_direct = A(i,k); % direct similarity
output.sim_common = p; % common neighbors score
output.sim_pref = (p + q)*(p + r); % preferential attachment index
output.sim_cos = p / sqrt( (p + q)*(p + r) ); % cosine coefficient
output.sim_Jaccard = p / (p + q + r); % Jaccard index
output.sim_Dice = 2*p / (2*p + q + r); % Dice coefficient
output.sim_prohub = p / min(p + q, p + r); % hub promoted index
output.sim_dehub = p / max(p + q, p + r); % hub depressed index

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
