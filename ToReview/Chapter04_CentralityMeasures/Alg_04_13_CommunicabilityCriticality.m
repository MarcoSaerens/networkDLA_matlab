function cr = Alg_04_13_CommunicabilityCriticality(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Masashi Shimbo (2018).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the communicability criticality of nodes.
%
% INPUT:
% ------- 
% - A: the (n x n) affinity matrix of an undirected graph.
%
% OUTPUT:
% -------
% - cr: (n x 1) vector of the communicability criticality scores.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n, nc] = size(A);
if n ~= nc
    error('affinity matrix not square');
elseif ~isequal(A, A')
    error('affinity matrix not symmetric');
end

K_c = expm(A);
K_c(K_c == 0) = 1; % replace 0 by 1 (necessary for non-connected graphs)

cr = zeros(n, 1);
for j = 1:n
    A_sans_j = A; 
    A_sans_j(j, :) = 0; A_sans_j(:, j) = 0; % set row & column j to 0
    
    K_c_sans_j = expm(A_sans_j); % communicability when node j is removed

    R_sans_j = 1 - K_c_sans_j ./ K_c;

    e_sans_j = ones(n, 1); e_sans_j(j) = 0; % e_sans_j = e - e_j
    cr(j) = e_sans_j' * (R_sans_j - diag(diag(R_sans_j))) * e_sans_j;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
