function ec = Alg_04_04_RandomEccentricity(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Masashi Shimbo
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the vectors of various centrality measures
%
% INPUT:
% ------- 
% - A : the (n x n) affinity matrix of a connected undirected graph
%
% OUTPUT:
% -------
% - ec : (n x 1) vector of random eccentricity scores of nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n, nc] = size(A);
if n ~= nc
  error('affinity matrix not square');
elseif A ~= A'
  error('affinity matrix not symmetric');
end

L = diag(sum(A, 1)) - A; % graph Laplacian

% P = inv(L + ones(n, n) / n) - ones(n, n) / n; % its pseudoinverse
P = pinv(L); % MATLAB's builtin pinv() is more stable

ec = diag(P); % eccentricity

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
