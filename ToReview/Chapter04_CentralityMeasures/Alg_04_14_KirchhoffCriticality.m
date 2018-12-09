function cr = Alg_04_14_KirchhoffCriticality(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Masashi Shimbo (2018).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Algorithm 4.14 for computing the node criticality based
%              on the Kirchhoff index.
%
% INPUT:
% ------- 
% - A: the (n x n) affinity matrix of an undirected connected graph.
%
% OUTPUT:
% -------
% - cr: (n x 1) vector of node criticality scores.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n, nc] = size(A);
if n ~= nc
    error('affinity matrix not square');
elseif ~isequal(A,  A')
    error('affinity matrix not symmetric');
end
% No check for connectedness; be careful

L = diag(sum(A, 1)) - A; % graph Laplacian
P = pinv(L); % pseudoinverse of L (= L^+)
P2 = P * P; % square of the pseudoinverse (= (L^+)^2)

cr = zeros(n, 1);
for j = 1:n
    % An implementation which is slow but more verbatim to the
    % pseudocode in the book would be:
    % for i = 1:n
    %     e_ij = zeros(n, 1); e_ij(i) = 1; e_ij(j) = -1;
    %     cr(j) = cr(j) + A(i,j) * e_ij' * P2 * e_ij;
    % end

    cr(j) = A(j,:) * ((diag(P2) + P2(j,j)) - 2 * P2(:,j));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
