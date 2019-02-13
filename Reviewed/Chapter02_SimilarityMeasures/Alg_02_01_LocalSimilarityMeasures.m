function sim = Alg_02_01_LocalSimilarityMeasures(A, i, k)

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
% - i, k : node indices to compare.
%
% OUTPUT:
% -------
% - sim : a structure containing various scalar similarities:
%              sim.direct  - direct similarity.
%              sim.common  - common neighbors score.
%              sim.pref    - preferential attachment index.
%              sim.cos     - cosine coefficient.
%              sim.Jaccard - Jaccard index.
%              sim.Dice    - Dice coefficient.
%              sim.prohub  - hub promoted index.
%              sim.dehub   - hub depressed index.
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

p = sum( A(i,:).*A(k,:));
q = sum( A(i,:).*(1 - A(k,:)) );
r = sum( (1 - A(i,:)).*A(k,:) );

sim.direct = A(i,k); % direct similarity
sim.common = p; % common neighbors score
sim.pref = (p + q)*(p + r); % preferential attachment index
sim.cos = p / sqrt( (p + q)*(p + r) ); % cosine coefficient
sim.Jaccard = p / (p + q + r); % Jaccard index
sim.Dice = 2*p / (2*p + q + r); % Dice coefficient
sim.prohub = p / min(p + q, p + r); % hub promoted index
sim.dehub = p / max(p + q, p + r); % hub depressed index

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
