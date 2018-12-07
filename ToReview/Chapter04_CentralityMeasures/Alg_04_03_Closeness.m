function [cc_k, cc_tc, cc_ks, cc_es] = Alg_04_03_Closeness(A, alpha_k, alpha_e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Various closeness centality measures
%
% Author: Masashi Shimbo
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Algorithm 4.3 of the above book for computing the vectors of
%              various centrality measures
%
% INPUT:
% ------- 
% - A : the (n x n) adjacency matrix of an undirected graph
% - alpha_k : discounting parameter (a non-negative real) for cc_k and cc_ks
%           : This value must be less than 1 / (spectral radius of A);
%           : otherwise, cc_k and cc_ks can be negative
% - alpha_e : discounting parameter for cc_tc and cc_es
%           : Unlike alpha_k, alpha_e can be any non-negative real number
%
% OUTPUT:
% -------
% - cc_k  : (n x 1) vector of the Katz centrality scores
% - cc_tc : (n x 1) vector of the total communicability centrality
% - cc_ks : (n x 1) vector of the Katz subgraph centrality
% - cc_es : (n x 1) vector of the exponential subgraph centrality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n, nc] = size(A);
if n ~= nc
  error('affinity matrix not square');
elseif A ~= A'
  error('adjacency matrix not symmetric');
end

%% Algorithm
e = ones(n, 1);
I = eye(n);
K = inv(I - alpha_k * A) - I; % Katz similarity
M = expm(alpha_e * A); % exponential diffusion kernel

% Katz centrality
cc_k = K * e;
% exponential (total communicability) centrality
cc_tc = M * e;
% Katz subgraph centrality
cc_ks = diag(K);
% exponential subgraph centrality
cc_es = diag(M);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
